"""Shared rendering of coincidence markers, used by BOTH the interactive app and
the plotting drivers so the on-page appearance is identical.

A coincidence is flagged with a small triangle just inside the panel edge
nearest the seam between the stacked EW/flux panels, using a blended transform
(x in data coords, y in axes fraction) so it stays put regardless of data
limits.  The marker *shape* encodes which row it is in; the *color* encodes the
kind:

  top row (EW)      -> '^' at y=edge          (10% above the bottom axis)
  bottom row (flux) -> 'v' at y=1-edge        (10% below the top axis)
  dip               -> C0 (blue),   alpha 0.3
  flare             -> C1 (orange), alpha 0.5

Call it once per axis, passing row='top' for the EW panel and row='bottom' for
the flux panel of a column.
"""
from __future__ import annotations


def draw_coincidence_markers(ax, dips=(), flares=(), *, row,
                             dip_color="C0", dip_alpha=0.3,
                             flare_color="C1", flare_alpha=0.5,
                             ms=4.0, edge=0.10, zorder=5):
    """Draw triangle markers at `dips` and `flares` (data x-coords) on `ax`.

    row='top'    -> '^' markers near the bottom of the panel
    row='bottom' -> 'v' markers near the top of the panel

    Returns the list of Line2D artists created (handy for the interactive app).
    """
    if row == "top":
        marker, y = "^", edge
    elif row == "bottom":
        marker, y = "v", 1.0 - edge
    else:
        raise ValueError(f"row must be 'top' or 'bottom', got {row!r}")

    trans = ax.get_xaxis_transform()  # x: data, y: axes fraction
    artists = []
    for xs, color, alpha in ((dips, dip_color, dip_alpha),
                             (flares, flare_color, flare_alpha)):
        xs = list(xs)
        if xs:
            (h,) = ax.plot(
                xs, [y] * len(xs), marker=marker, linestyle="none",
                color=color, alpha=alpha, markersize=ms, markeredgewidth=0,
                transform=trans, clip_on=False, zorder=zorder,
            )
            artists.append(h)
    return artists
