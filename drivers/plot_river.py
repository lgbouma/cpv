"""
Given a TIC ID, make a river plot.
"""

import os
import numpy as np, pandas as pd
from complexrotators.getters import _get_lcpaths_fromlightkurve_given_ticid

from complexrotators.plotting import plot_river
from complexrotators.paths import CSVDIR, RESULTSDIR, LOCALDIR
from os.path import join

from complexrotators.lcprocessing import (
    prepare_cpv_light_curve
)


cmap = 'seismic'
#cmap = 'Paired'
#cmap = 'gist_stern'
#cmap = 'PiYG'

def many_ticids():

    raise NotImplementedError('deprecated')

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

    ##########
    # TIC 4029 CONFIGURATION
    ticid = '402980664'
    cyclewindows = [(0,1481), (0,64), (248,315), (1233, 1265), (1233,1481), (1411, 1481)]
    sample_id = '2023catalog_LGB_RJ_concat' # used for cacheing
    manual_period = 18.5611 / 24
    t0 = 1791.15 + manual_period/2  # can be None
    # e.g. written by build_4029_mask.py
    nterms = 2
    model_id = f"manual_20230617_mask_v0_nterms{nterms}"
    # specifies lc subtrxn model
    manual_csvpath = f'/Users/luke/Dropbox/proj/cpv/results/4029_mask/lc_lsresid_{model_id}.csv'

    savstr = model_id # can be None; used for cacheing

    vmin, vmax = -0.03, 0.03 #  can be None
    vmin, vmax = -0.04, 0.01 #  can be None
    ##########

    ##########
    # TIC 405910546
    ticid = '405910546'
    cyclewindows = [(-11,27)] #[(0,1481), (0,64), (248,315), (1233, 1265), (1233,1481), (1411, 1481)]
    sample_id = '2023catalog_LGB_RJ_concat' # used for cacheing
    manual_period = 1.583
    t0 = 2321.77 + manual_period/2
    # e.g. written by build_4029_mask.py
    nterms = None
    model_id = None
    # specifies lc subtrxn model
    manual_csvpath = None
    savstr = model_id # can be None; used for cacheing
    vmin, vmax = -0.07, 0.05 #  can be None
    ##########

    ##########
    ticid = '300651846'
    cyclewindows = [(0, 650), (2305, 2630)]
    sample_id = '2023catalog_LGB_RJ_concat'
    manual_period =  8.2551/24
    manual_period =  8.254/24  #better
    t0 = 2170. - 5*manual_period
    nterms = None
    model_id = None
    # specifies lc subtrxn model
    manual_csvpath = None
    savstr = model_id # can be None; used for cacheing
    vmin, vmax = -0.06, 0.06 #  can be None

    ##########################################

    # directory management
    outdir = join(RESULTSDIR, 'river')
    if not os.path.exists(outdir): os.mkdir(outdir)
    outdir = join(RESULTSDIR, 'river', f'TIC_{ticid}')
    if not os.path.exists(outdir): os.mkdir(outdir)

    cachedir = join(LOCALDIR, "cpv_finding")
    if not os.path.exists(cachedir): os.mkdir(cachedir)
    cachedir = join(cachedir, sample_id)
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    if isinstance(manual_csvpath, str):
        df = pd.read_csv(manual_csvpath)
        times = np.array(df.time)
        # residual flux from subtracting the model given in model_id
        fluxs = np.array(df.r_flux)

    else:
        lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid)

        # stack over all sectors
        times, fluxs = [], []
        for lcpath in np.sort(lcpaths):
            (_, _, _, x_obs, _, y_flat, _, _, _,
             _, _) = prepare_cpv_light_curve(lcpath, cachedir)
            times.append(x_obs)
            fluxs.append(y_flat)
        times = np.hstack(np.array(times).flatten())
        fluxs = np.hstack(np.array(fluxs).flatten())


    idstr = f'TIC{ticid}'
    titlestr = f'TIC{ticid}. P = {manual_period*24:.4f}hr.'

    # make the plots
    for cyclewindow in cyclewindows:

        plot_river(times, fluxs, manual_period, outdir, titlestr=titlestr,
                   cmap=cmap, cyclewindow=cyclewindow, idstr=idstr, t0=t0,
                   savstr=savstr, vmin=vmin, vmax=vmax)



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
