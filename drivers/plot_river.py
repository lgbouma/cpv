"""
Given a TIC ID, make a river plot.
"""

from complexrotators.helpers import (
    get_complexrot_data,
    get_complexrot_twentysec_data
)
from complexrotators.plotting import plot_river
from complexrotators.paths import CSVDIR
from os.path import join

cmap = 'seismic'

def many_ticids():

    #ticids = ['177309964','206544316','300651846','201789285']
    #cyclewindowss = [
    #    [None, (0,800), (1550,2000)],
    #    [None, (0,180), (2280, 2450)],
    #    [None],
    #    [None,(0,350),(4830,5200)]
    #]
    ticids = '59129133,245902096,118769116'.split(',')
    cyclewindowss = [[None],[None],[None]]

    for ticid, cyclewindows in zip(ticids, cyclewindowss):

        d = get_complexrot_data(ticid)

        idstr = f'TIC{ticid}'
        titlestr = f'TIC{ticid}. P: {d["period"]:.7f}d'

        for cyclewindow in cyclewindows:
            plot_river(d['times'], d['fluxs'], d['period'], d['outdir'],
                       titlestr=titlestr, cmap=cmap, cyclewindow=cyclewindow,
                       idstr=idstr)


def single_ticid():

    ##########################################
    # change these
    ticid = '364075855'
    cyclewindows = [None,(0,65),(145,190)] # [None] #
    is_20s = False
    MANUAL_PERIOD = 1.317741

    # typically None
    readcsv = join(
        CSVDIR,
        'TIC364075855_cdips_tess_lightcurve_minimal_PCA1_allsector_sigclipped.csv'
    )
    ##########################################

    # important keys: times, fluxs, period, t0, lsp.
    if is_20s:
        d = get_complexrot_twentysec_data(ticid)
    else:
        d = get_complexrot_data(ticid, hardcsv=readcsv)

    d['period'] = MANUAL_PERIOD

    idstr = f'TIC{ticid}'
    titlestr = f'TIC{ticid}. P: {d["period"]:.7f}d'

    for cyclewindow in cyclewindows:

        plot_river(d['times'], d['fluxs'], d['period'], d['outdir'],
                   titlestr=titlestr, cmap=cmap, cyclewindow=cyclewindow,
                   idstr=idstr)



def single_kicid():

    # kicid = '7740983' # saul's Kepler CR
    kicid = '6184894' # kepler1627
    cyclewindows = [None] #[None,(0,350),(4830,5200)]

    # important keys: times, fluxs, period, t0, lsp.
    d = get_complexrot_data(None, kicid=kicid)

    idstr = f'KIC{kicid}'
    titlestr = f'KIC{kicid}. P: {d["period"]:.7f}d'

    for cyclewindow in cyclewindows:

        plot_river(d['times'], d['fluxs'], d['period'], d['outdir'],
                   titlestr=titlestr, cmap=cmap, cyclewindow=cyclewindow,
                   idstr=idstr)




if __name__ == "__main__":

    single_ticid()
    # single_kicid()
    #many_ticids()
