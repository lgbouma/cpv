"""Interactive matplotlib labeler for one mosaic figure.

Reproduces the figure's columns (EW on top, Δflux on bottom, shared x) and lets
you place / remove coincidence markers by clicking:

  left-click   in a column -> add a DIP  ('|')  at the nearest EW datapoint x
  right-click  in a column -> add a FLARE ('v') at the nearest EW datapoint x
  click on/near an existing marker (either button) -> remove it
  'u'          undo the last add/remove
  'w'          write labels now (also auto-saved after every change)
  'q'          close this figure / advance

Markers are drawn above the top (EW) panel; because the panels share an x-axis,
they flag the aligned location in the flux panel too.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .data import figure_columns, build_column_series
from .labels import (colkey, load_labels, save_labels, get_column_labels,
                     _COLUMNS)
from .render import draw_coincidence_markers

_TOL_PX = 9.0  # click-to-marker match tolerance, pixels (x only)


class LabelSession:
    def __init__(self, figure: str):
        self.figure = figure
        self.cols = figure_columns(figure)
        self.ncols = len(self.cols)
        self.series = [build_column_series(figure, c) for c in self.cols]
        self.colkeys = [colkey(c) for c in self.cols]

        # per-column state: {'dip': [x...], 'flare': [x...]}
        df = load_labels()
        self.state = []
        for ck in self.colkeys:
            self.state.append({
                "dip": list(get_column_labels(df, figure, ck, "dip")),
                "flare": list(get_column_labels(df, figure, ck, "flare")),
            })

        self.history = []            # (col_ix, kind, x, action) for undo
        self.marker_artists = [[] for _ in range(self.ncols)]
        self.ax_to_col = {}
        self.ew_axes = []

        self._build_figure()

    # ── figure construction ────────────────────────────────────────────────
    def _build_figure(self):
        _F = 0.9
        self.fig, self.axs = plt.subplots(
            nrows=2, ncols=self.ncols,
            figsize=(_F * self.ncols * 2.4, _F * 4.6),
            sharex="col", squeeze=False,
        )
        c = "k"
        for i in range(self.ncols):
            ax_ew = self.axs[0, i]
            ax_fl = self.axs[1, i]
            self.ax_to_col[ax_ew] = i
            self.ax_to_col[ax_fl] = i
            self.ew_axes.append(ax_ew)

            s = self.series[i]
            ax_ew.set_title(s["label"], fontsize="x-small", pad=3)
            ax_ew.errorbar(
                s["ew"]["x"], s["ew"]["y"], yerr=s["ew"]["yerr"],
                fmt="o", color=c, ecolor=c, elinewidth=0.5,
                capsize=0, markersize=1.8,
            )
            ax_fl.errorbar(
                s["flux"]["x"], s["flux"]["y"], yerr=s["flux"]["yerr"],
                fmt="o", color=c, ecolor=c, elinewidth=0.2,
                capsize=0, markersize=1.2, mew=0,
            )
            ax_fl.set_ylim(*s["flux_ylim"])
            if s.get("ew_ymin") is not None:
                ax_ew.set_ylim(bottom=s["ew_ymin"])
            ax_ew.set_xlim(*s["xlim"])
            ax_fl.set_xlabel(s["xlabel"], fontsize="small")

        self.axs[0, 0].set_ylabel(r"H$\alpha$ EW [$\mathrm{\AA}$]",
                                  fontsize="small")
        self.axs[1, 0].set_ylabel(r"$\Delta$ Flux [%]", fontsize="small")
        for ax in self.axs.flatten():
            ax.tick_params(axis="both", which="major", labelsize="x-small")

        self.fig.suptitle(
            f"[{self.figure}]  L-click: dip (blue)   R-click: flare (orange)"
            "   click marker: remove   u: undo   w: save   q: done",
            fontsize="small",
        )
        self.fig.tight_layout(rect=(0, 0, 1, 0.96))

        for i in range(self.ncols):
            self._redraw_col(i)

        self.fig.canvas.mpl_connect("button_press_event", self._on_click)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

    # ── drawing ─────────────────────────────────────────────────────────────
    def _redraw_col(self, i):
        for h in self.marker_artists[i]:
            h.remove()
        dips, flares = self.state[i]["dip"], self.state[i]["flare"]
        arts = []
        arts += draw_coincidence_markers(self.axs[0, i], dips, flares, row="top")
        arts += draw_coincidence_markers(self.axs[1, i], dips, flares, row="bottom")
        self.marker_artists[i] = arts
        self.fig.canvas.draw_idle()

    # ── event handling ──────────────────────────────────────────────────────
    def _toolbar_active(self):
        tb = getattr(self.fig.canvas, "toolbar", None)
        return bool(tb and getattr(tb, "mode", ""))

    def _nearest_existing(self, i, ax, x):
        """(kind, x, px_dist) of the closest existing marker in column i, or None."""
        px_click = ax.transData.transform((x, 0))[0]
        best = None
        for kind in ("dip", "flare"):
            for xv in self.state[i][kind]:
                px = ax.transData.transform((xv, 0))[0]
                d = abs(px - px_click)
                if best is None or d < best[2]:
                    best = (kind, xv, d)
        return best

    def _snap(self, i, x):
        snap = self.series[i]["snap_x"]
        return float(snap[np.argmin(np.abs(snap - x))])

    def _on_click(self, event):
        if self._toolbar_active():
            return
        if event.inaxes is None or event.inaxes not in self.ax_to_col:
            return
        if event.xdata is None:
            return
        i = self.ax_to_col[event.inaxes]
        ax = event.inaxes

        # Remove if the click lands on an existing marker.
        near = self._nearest_existing(i, ax, event.xdata)
        if near is not None and near[2] <= _TOL_PX:
            kind, xv, _ = near
            self.state[i][kind].remove(xv)
            self.history.append((i, kind, xv, "remove"))
            self._redraw_col(i)
            self._save()
            return

        # Otherwise add: left -> dip, right -> flare.
        if event.button == 1:
            kind = "dip"
        elif event.button == 3:
            kind = "flare"
        else:
            return
        xv = self._snap(i, event.xdata)
        if xv in self.state[i][kind]:
            return
        self.state[i][kind].append(xv)
        self.history.append((i, kind, xv, "add"))
        self._redraw_col(i)
        self._save()

    def _on_key(self, event):
        if event.key == "u":
            self._undo()
        elif event.key == "w":
            self._save()
            print(f"[{self.figure}] labels saved.")
        elif event.key == "q":
            plt.close(self.fig)

    def _undo(self):
        if not self.history:
            return
        i, kind, xv, action = self.history.pop()
        if action == "add":
            if xv in self.state[i][kind]:
                self.state[i][kind].remove(xv)
        else:  # remove -> re-add
            self.state[i][kind].append(xv)
        self._redraw_col(i)
        self._save()

    # ── persistence ─────────────────────────────────────────────────────────
    def _save(self):
        df = load_labels()
        df = df[df.figure != self.figure]
        rows = []
        for i, ck in enumerate(self.colkeys):
            for kind in ("dip", "flare"):
                for x in self.state[i][kind]:
                    rows.append({"figure": self.figure, "colkey": ck,
                                 "kind": kind, "x": x})
        new = pd.DataFrame(rows, columns=_COLUMNS)
        save_labels(pd.concat([df, new], ignore_index=True))

    def run(self):
        plt.show()
