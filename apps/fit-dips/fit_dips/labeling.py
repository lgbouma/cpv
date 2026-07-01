"""
Interactive matplotlib labeling UI for CPV dip fitting.

LabelState is a headless, unit-testable container of the user's dip/flare
windows.  DipLabeler is the matplotlib GUI that drives the workflow:

    drag spans -> dip windows        (default mode)
    'f' then drag -> flare windows   ('d' switches back to dip mode)
    'u' undo last span, 'r' reset
    'enter' -> classify points + ask to confirm
    'y' -> run fits, overlay preferred model      'n' -> keep editing
    'q' -> quit

The fitting itself is injected via the ``on_fit`` callback so this module
stays decoupled from fitting/registry.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector


class LabelState:
    """Headless container for dip/flare time-windows."""

    def __init__(self):
        self.dips = []     # list of [x0, x1]
        self.flares = []   # list of [x0, x1]
        self.order = []    # ('dip'|'flare', index) for undo
        self.confirmed = False

    def add_dip(self, x0, x1):
        self.dips.append([float(min(x0, x1)), float(max(x0, x1))])
        self.order.append(("dip", len(self.dips) - 1))

    def add_flare(self, x0, x1):
        self.flares.append([float(min(x0, x1)), float(max(x0, x1))])
        self.order.append(("flare", len(self.flares) - 1))

    def undo(self):
        if not self.order:
            return
        kind, _ = self.order.pop()
        if kind == "dip":
            self.dips.pop()
        else:
            self.flares.pop()

    def reset(self):
        self.dips, self.flares, self.order = [], [], []
        self.confirmed = False

    def in_windows(self, t, windows):
        t = np.asarray(t, dtype=float)
        m = np.zeros(t.size, dtype=bool)
        for (x0, x1) in windows:
            m |= (t >= x0) & (t <= x1)
        return m

    def masks(self, t):
        """Return (in_dip, in_flare, out_of_dip) boolean masks."""
        in_dip = self.in_windows(t, self.dips)
        in_flare = self.in_windows(t, self.flares)
        in_dip = in_dip & ~in_flare          # flare wins over dip
        oot = ~in_dip & ~in_flare
        return in_dip, in_flare, oot


_INSTRUCTIONS = (
    "drag=mark dip | 'x'=flare mode | 'd'=dip mode | 'u'=undo | "
    "'r'=reset | enter=classify | 'y'=fit | 'n'=edit | 'q'=quit"
)


def _phase_bin(t, flux, period, per_period=200):
    """Bin (t, flux) into ~``per_period`` equal-width bins per ``period``.

    Intended for a phase-folded light curve (times already folded onto one
    cycle, e.g. the _S88 variant), so the data span is ~one period and the
    result is ~``per_period`` binned points. Returns (t_mean, flux_mean) over
    the non-empty bins. Falls back to the raw arrays if the span/period is
    degenerate.
    """
    t = np.asarray(t, dtype=float)
    flux = np.asarray(flux, dtype=float)
    lo, hi = float(t.min()), float(t.max())
    span = hi - lo
    if span <= 0 or not np.isfinite(period) or period <= 0:
        return t, flux
    nb = max(1, int(np.ceil(per_period * span / period)))
    edges = np.linspace(lo, hi, nb + 1)
    cnt, _ = np.histogram(t, bins=edges)
    tsum, _ = np.histogram(t, bins=edges, weights=t)
    fsum, _ = np.histogram(t, bins=edges, weights=flux)
    m = cnt > 0
    return tsum[m] / cnt[m], fsum[m] / cnt[m]


class MultiDipLabeler:
    """Stacked, shared-x labeler over one or more co-temporal light curves.

    Each panel is a dict ``{t, flux, title, on_fit, bin_period}``. All panels
    share a single :class:`LabelState`, so a dip/flare window dragged in any
    panel propagates to every panel -- guaranteeing identical in/out-of-dip
    times across bands on the same night. The masks are still evaluated per
    panel against that panel's own time array, so different per-band sampling
    is handled correctly. Pressing 'y' fits *every* panel through its own
    ``on_fit`` callback and overlays each preferred model.

    The single-panel case (``len(panels) == 1``) reproduces the original
    one-axis labeler, including the optional phase-binned overlay; see
    :class:`DipLabeler` for the thin one-panel wrapper.
    """

    def __init__(self, panels, group_title, on_group_fit=None):
        assert len(panels) >= 1
        self.group_title = group_title
        # Optional group-level fit hook: when set, 'y' calls it ONCE with the
        # shared LabelState and expects one overlay payload per panel (used for
        # the cross-band shared-duration joint fit). When None, each panel is
        # fit independently through its own on_fit (single-band / --solo).
        self.on_group_fit = on_group_fit
        self.state = LabelState()
        self.mode = "dip"
        self.awaiting_confirm = False
        self._spans = []   # axvspan patches across ALL axes

        n = len(panels)
        height = 6.0 if n == 1 else (2.7 * n + 1.4)
        self.fig, axes = plt.subplots(
            n, 1, sharex=True, figsize=(11, height), squeeze=False)
        self.axes = list(axes[:, 0])

        # Per-panel render state lives in the panel dict alongside its data.
        self.panels = []
        all_t = []
        for ax, p in zip(self.axes, panels):
            panel = dict(p)
            panel["ax"] = ax
            panel["t"] = np.asarray(panel["t"], dtype=float)
            panel["flux"] = np.asarray(panel["flux"], dtype=float)
            panel.setdefault("bin_period", None)
            panel["binned"] = None
            panel["class_artists"] = []   # colored points from classify preview
            panel["model_line"] = None
            panel["baseline_line"] = None
            panel["dip_marker_lines"] = []
            all_t.append(panel["t"])
            self._render_base(panel)
            ax.set_ylabel("normalized flux")
            if panel.get("title"):
                ax.set_title(panel["title"], fontsize=8, loc="left")
            self.panels.append(panel)
        self.axes[-1].set_xlabel("time [BJD]")

        # Shared, tight x-range over the union of all panels' times.
        tmin = min(float(t.min()) for t in all_t)
        tmax = max(float(t.max()) for t in all_t)
        span = tmax - tmin
        pad = 0.02 * span if span > 0 else 1.0
        self.axes[0].set_xlim(tmin - pad, tmax + pad)
        for ax in self.axes:
            ax.margins(x=0)  # keep span/overlay redraws from re-padding x

        # One SpanSelector per axis (kept referenced so they aren't GC'd); each
        # routes to the shared handler so a drag in any panel labels them all.
        self._selectors = [
            SpanSelector(
                ax, self._on_select, "horizontal", useblit=True,
                props=dict(alpha=0.25, facecolor="tab:blue"),
                interactive=False, drag_from_anywhere=True)
            for ax in self.axes]

        # Disconnect matplotlib's default key bindings so they don't shadow ours
        # (e.g. 'f'=fullscreen, 'r'/'h'=home, 's'=save, 'p'=pan, 'o'=zoom).
        # Toolbar buttons for pan/zoom/save/home still work.
        try:
            self.fig.canvas.mpl_disconnect(
                self.fig.canvas.manager.key_press_handler_id)
        except Exception:
            pass
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

        self.fig.subplots_adjust(top=0.84, hspace=0.18)
        self._suptitle = None
        self._set_title()

    # -- display helpers --------------------------------------------------
    def _render_base(self, panel):
        """Draw the raw (and optional phase-binned) points for one panel."""
        ax, t, flux = panel["ax"], panel["t"], panel["flux"]
        if panel["bin_period"]:
            # Dense folded light curve (e.g. _S88): faint gray raw points with
            # a black ~200/period phase-binned overlay, so dip/flare edges stay
            # legible while dragging out the regions to fit.
            panel["pts"] = ax.plot(t, flux, ".", color="0.6", ms=2,
                                   alpha=0.5, zorder=2)[0]
            tb, fb = _phase_bin(t, flux, panel["bin_period"])
            panel["binned"] = ax.plot(
                tb, fb, "o", color="k", ms=4, alpha=1.0, mec="white", mew=0.4,
                zorder=4, label="binned (200/period)")[0]
            ax.legend(loc="best", fontsize=8)
        else:
            panel["pts"] = ax.plot(t, flux, ".", color="k", ms=3, zorder=3)[0]

    def _set_title(self, extra=""):
        mode_txt = f"MODE={self.mode.upper()}"
        ndip = len(self.state.dips)
        nfl = len(self.state.flares)
        self._suptitle = self.fig.suptitle(
            f"{self.group_title}\n{mode_txt}  dips={ndip} flares={nfl}  "
            f"{extra}\n{_INSTRUCTIONS}", fontsize=9)
        self.fig.canvas.draw_idle()

    def _draw_span(self, x0, x1, color):
        for ax in self.axes:
            self._spans.append(ax.axvspan(x0, x1, color=color, alpha=0.18,
                                          zorder=0))

    def _redraw_spans(self):
        for s in self._spans:
            s.remove()
        self._spans = []
        for (x0, x1) in self.state.dips:
            self._draw_span(x0, x1, "tab:blue")
        for (x0, x1) in self.state.flares:
            self._draw_span(x0, x1, "tab:red")

    # -- event handlers ---------------------------------------------------
    def _on_select(self, x0, x1):
        if x1 - x0 <= 0:
            return
        if self.mode == "dip":
            self.state.add_dip(x0, x1)
            self._draw_span(x0, x1, "tab:blue")
        else:
            self.state.add_flare(x0, x1)
            self._draw_span(x0, x1, "tab:red")
        self.awaiting_confirm = False
        self._set_title()

    def _classify_preview(self):
        n_oot_total = 0
        for panel in self.panels:
            ax, t, flux = panel["ax"], panel["t"], panel["flux"]
            in_dip, in_flare, oot = self.state.masks(t)
            n_oot_total += int(np.sum(oot))
            # Clear any prior preview (raw points + last classification).
            if panel["pts"] is not None:
                panel["pts"].remove()
                panel["pts"] = None
            for a in panel["class_artists"]:
                a.remove()
            panel["class_artists"] = []
            panel["class_artists"] += [
                ax.plot(t[oot], flux[oot], ".", color="k", ms=3, zorder=3,
                        label="out-of-dip")[0],
                ax.plot(t[in_dip], flux[in_dip], ".", color="tab:blue", ms=4,
                        zorder=4, label="in-dip")[0],
                ax.plot(t[in_flare], flux[in_flare], ".", color="tab:red",
                        ms=4, zorder=4, label="flare")[0]]
            ax.legend(loc="best", fontsize=8)
        self.awaiting_confirm = True
        self._set_title(
            extra=f"CONFIRM? {n_oot_total} pts auto out-of-dip -> 'y'/'n'")

    def _overlay_fit(self, panel, payload):
        ax = panel["ax"]
        if panel["model_line"] is not None:
            panel["model_line"].remove()
        if panel["baseline_line"] is not None:
            panel["baseline_line"].remove()
            panel["baseline_line"] = None
        for ln in panel["dip_marker_lines"]:
            ln.remove()
        panel["dip_marker_lines"] = []
        panel["model_line"] = ax.plot(
            payload["tgrid"], payload["model"], "-", color="tab:green", lw=2,
            zorder=5, label=f"preferred: {payload['preferred']}")[0]
        if "baseline" in payload:
            panel["baseline_line"] = ax.plot(
                payload["tgrid"], payload["baseline"], ":", color="tab:orange",
                lw=1.5, zorder=5, label="baseline (no dip)")[0]
        for t0 in payload.get("dip_markers", []):
            panel["dip_marker_lines"].append(
                ax.axvline(t0, color="tab:green", ls=":", lw=1, zorder=2))
        ax.legend(loc="best", fontsize=8)

    def _do_fit(self):
        self._set_title(extra="fitting all panels...")
        plt.pause(0.01)
        if self.on_group_fit is not None:
            # One joint fit over the whole group; payloads align with panels.
            payloads = self.on_group_fit(self.state)
        else:
            payloads = [panel["on_fit"](self.state) for panel in self.panels]
        preferred = []
        for panel, payload in zip(self.panels, payloads):
            if payload is None:
                continue
            self._overlay_fit(panel, payload)
            preferred.append(payload["preferred"])
            header = panel.get("title") or "panel"
            print(f"\n=== {header} ===")
            print(payload["table_text"])
        self.awaiting_confirm = False
        self.state.confirmed = True
        n = len(preferred)
        self._set_title(extra=f"DONE -> fit {n} panel(s) (saved). "
                              f"'r' to redo, 'q' to quit.")

    def _clear_fit_overlays(self):
        for panel in self.panels:
            if panel["model_line"] is not None:
                panel["model_line"].remove()
                panel["model_line"] = None
            if panel["baseline_line"] is not None:
                panel["baseline_line"].remove()
                panel["baseline_line"] = None
            for ln in panel["dip_marker_lines"]:
                ln.remove()
            panel["dip_marker_lines"] = []

    def _on_key(self, event):
        k = event.key
        if k == "x":
            self.mode = "flare"
            for s in self._selectors:
                s.set_props(facecolor="tab:red")
            self._set_title()
        elif k == "d":
            self.mode = "dip"
            for s in self._selectors:
                s.set_props(facecolor="tab:blue")
            self._set_title()
        elif k == "u":
            self.state.undo()
            self._redraw_spans()
            self._set_title()
        elif k == "r":
            self.state.reset()
            self._redraw_spans()
            self._clear_fit_overlays()
            self._set_title(extra="reset")
        elif k == "enter":
            if len(self.state.dips) == 0:
                self._set_title(extra="mark at least one dip first")
                return
            self._classify_preview()
        elif k == "y":
            if self.awaiting_confirm:
                self._do_fit()
        elif k == "n":
            if self.awaiting_confirm:
                self.awaiting_confirm = False
                self._set_title(extra="keep editing")
        elif k == "q":
            plt.close(self.fig)

    def run(self):
        plt.show()


class DipLabeler(MultiDipLabeler):
    """Single-panel labeler -- a thin one-panel :class:`MultiDipLabeler`."""

    def __init__(self, t, flux, title, on_fit, bin_period=None):
        super().__init__(
            [{"t": t, "flux": flux, "title": "", "on_fit": on_fit,
              "bin_period": bin_period}],
            group_title=title)
        # Back-compat aliases for the original single-axis attribute names.
        self.ax = self.axes[0]
