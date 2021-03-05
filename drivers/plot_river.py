"""
Given a TIC ID, make a river plot.
"""

from complexrotators.helpers import get_complexrot_data
from complexrotators.plotting import plot_river

ticids = ['177309964','206544316','300651846','201789285']
cyclewindowss = [
    [None, (0,800), (1550,2000)],
    [None, (0,180), (2280, 2450)],
    [None],
    [None,(0,350),(4830,5200)]
]

def many_ticids():

    # cmap = 'viridis'
    # cmap = 'bwr'
    # cmap = 'BrBG'
    cmap = 'spectral'
    cmap = 'seismic'

    for ticid, cyclewindows in zip(ticids, cyclewindowss):

        d = get_complexrot_data(ticid)

        idstr = f'TIC{ticid}'
        titlestr = f'TIC{ticid}. P: {d["period"]:.7f}d'

        for cyclewindow in cyclewindows:
            plot_river(d['times'], d['fluxs'], d['period'], d['outdir'],
                       titlestr=titlestr, cmap=cmap, cyclewindow=cyclewindow,
                       idstr=idstr)


def single_ticid():

    ticid = '201789285'
    cyclewindows = [None,(0,350),(4830,5200)]

    # important keys: times, fluxs, period, t0, lsp.
    d = get_complexrot_data(ticid)

    idstr = f'TIC{ticid}'
    titlestr = f'TIC{ticid}. P: {d["period"]:.7f}d'

    for cyclewindow in cyclewindows:

        plot_river(d['times'], d['fluxs'], d['period'], d['outdir'],
                   titlestr=titlestr, cmap='viridis', cyclewindow=cyclewindow,
                   idstr=idstr)


if __name__ == "__main__":
    many_ticids()
