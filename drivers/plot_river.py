"""
Given a TIC ID, make a river plot.
"""

from complexrotators.helpers import get_complexrot_data
from complexrotators.plotting import plot_river

ticid = '177309964'
cyclewindows = [None, (0,800), (1550,2000)]

# important keys: times, fluxs, period, t0, lsp.
d = get_complexrot_data(ticid)

titlestr = f'TIC{ticid}'

for cyclewindow in cyclewindows:

    plot_river(d['times'], d['fluxs'], d['period'], d['outdir'],
               titlestr=titlestr, cmap='viridis', cyclewindow=cyclewindow)
