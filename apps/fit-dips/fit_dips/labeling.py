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


class DipLabeler:
    def __init__(self, t, flux, title, on_fit):
        self.t = np.asarray(t, dtype=float)
        self.flux = np.asarray(flux, dtype=float)
        self.title = title
        self.on_fit = on_fit
        self.state = LabelState()
        self.mode = "dip"
        self.awaiting_confirm = False
        self._spans = []   # axvspan patches

        self.fig, self.ax = plt.subplots(figsize=(11, 6))
        self.pts = self.ax.plot(self.t, self.flux, ".", color="k", ms=3,
                                zorder=3)[0]
        self.model_line = None
        self.baseline_line = None
        # Tight default x-range: fill the axes with the data plus a small pad.
        span = float(self.t.max() - self.t.min())
        pad = 0.02 * span if span > 0 else 1.0
        self.ax.set_xlim(self.t.min() - pad, self.t.max() + pad)
        self.ax.margins(x=0)  # keep span/overlay redraws from re-padding x
        self.ax.set_xlabel("time [BJD]")
        self.ax.set_ylabel("normalized flux")
        self._set_title()

        self.span = SpanSelector(
            self.ax, self._on_select, "horizontal", useblit=True,
            props=dict(alpha=0.25, facecolor="tab:blue"),
            interactive=False, drag_from_anywhere=True)
        # Disconnect matplotlib's default key bindings so they don't shadow ours
        # (e.g. 'f'=fullscreen, 'r'/'h'=home, 's'=save, 'p'=pan, 'o'=zoom).
        # Toolbar buttons for pan/zoom/save/home still work.
        try:
            self.fig.canvas.mpl_disconnect(
                self.fig.canvas.manager.key_press_handler_id)
        except Exception:
            pass
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

    # -- display helpers --------------------------------------------------
    def _set_title(self, extra=""):
        mode_txt = f"MODE={self.mode.upper()}"
        ndip = len(self.state.dips)
        nfl = len(self.state.flares)
        self.ax.set_title(
            f"{self.title}\n{mode_txt}  dips={ndip} flares={nfl}  {extra}\n"
            f"{_INSTRUCTIONS}", fontsize=9)
        self.fig.canvas.draw_idle()

    def _draw_span(self, x0, x1, color):
        self._spans.append(
            self.ax.axvspan(x0, x1, color=color, alpha=0.18, zorder=0))

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
        in_dip, in_flare, oot = self.state.masks(self.t)
        self.pts.remove()
        self.ax.plot(self.t[oot], self.flux[oot], ".", color="k", ms=3,
                     zorder=3, label="out-of-dip")
        self.ax.plot(self.t[in_dip], self.flux[in_dip], ".", color="tab:blue",
                     ms=4, zorder=4, label="in-dip")
        self.ax.plot(self.t[in_flare], self.flux[in_flare], ".",
                     color="tab:red", ms=4, zorder=4, label="flare")
        self.pts = self.ax.plot([], [], ".", color="k", ms=3)[0]
        self.ax.legend(loc="best", fontsize=8)
        self.awaiting_confirm = True
        n_oot = int(np.sum(oot))
        self._set_title(extra=f"CONFIRM? {n_oot} pts auto out-of-dip -> 'y'/'n'")

    def _do_fit(self):
        self._set_title(extra="fitting...")
        plt.pause(0.01)
        payload = self.on_fit(self.state)
        if payload is None:
            self._set_title(extra="fit returned nothing")
            return
        tgrid = payload["tgrid"]
        model = payload["model"]
        if self.model_line is not None:
            self.model_line.remove()
        if self.baseline_line is not None:
            self.baseline_line.remove()
            self.baseline_line = None
        self.model_line = self.ax.plot(
            tgrid, model, "-", color="tab:green", lw=2, zorder=5,
            label=f"preferred: {payload['preferred']}")[0]
        if "baseline" in payload:
            self.baseline_line = self.ax.plot(
                tgrid, payload["baseline"], ":", color="tab:orange", lw=1.5,
                zorder=5, label="baseline (no dip)")[0]
        for t0 in payload.get("dip_markers", []):
            self.ax.axvline(t0, color="tab:green", ls=":", lw=1, zorder=2)
        self.ax.legend(loc="best", fontsize=8)
        self.awaiting_confirm = False
        self.state.confirmed = True
        print("\n" + payload["table_text"])
        self._set_title(extra=f"DONE -> preferred {payload['preferred']} "
                              f"(saved). 'r' to redo, 'q' to quit.")

    def _on_key(self, event):
        k = event.key
        if k == "x":
            self.mode = "flare"
            self.span.set_props(facecolor="tab:red")
            self._set_title()
        elif k == "d":
            self.mode = "dip"
            self.span.set_props(facecolor="tab:blue")
            self._set_title()
        elif k == "u":
            self.state.undo()
            self._redraw_spans()
            self._set_title()
        elif k == "r":
            self.state.reset()
            self._redraw_spans()
            if self.model_line is not None:
                self.model_line.remove()
                self.model_line = None
            if self.baseline_line is not None:
                self.baseline_line.remove()
                self.baseline_line = None
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
