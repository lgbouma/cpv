"""
This is a simple pipeline for identifying CPVs.

Contents:
    | find_CPV: given ticid, determine if it's a CPV
    | find_CPVs: thin wrapper
    | prepare_cpv_light_curve: retrieve all relevant data from SPOC/FITS LC
"""
import os, pickle, subprocess
from os.path import join
import numpy as np, pandas as pd
from glob import glob

from complexrotators.paths import LOCALDIR
from complexrotators.plotting import (
    plot_phased_light_curve, plot_dipcountercheck
)
from complexrotators.getters import _get_lcpaths_given_ticid
from complexrotators.lcprocessing import (
    cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
)

def find_CPVs():

    # the TICIDs to search
    ticids = [
    "201789285",
    "311092148",
    "332517282",
    "405910546",
    "142173958",
    "300651846",
    "408188366",
    "146539195",
    "177309964",
    "425933644",
    "206544316",
    "224283342",
    "245902096",
    "150068381",
    "177309964",
    "118769116",
    "245868207",
    "245874053",
    "59129133"
    ]

    for ticid in ticids:
        find_CPV(ticid)


def find_CPV(ticid):

    #
    # get the light curves for all desired sectors and cadences
    #
    lcpaths = _get_lcpaths_given_ticid(ticid)

    cachedir = join(LOCALDIR, "cpv_finding")
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    #
    # for each light curve (sector / cadence specific), detrend, ((remove
    # flares)), get the best period, and then phase-fold.
    #
    for lcpath in lcpaths:

        logpath = join(cachedir, f'{starid}.log')

        # get the relevant light curve data
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        # get period, t0, and periodogram (PDM or LombScargle)
        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        # require PDM theta statistic < 0.8
        # and period < 2 days
        #TODO TODO
        import IPython; IPython.embed()
        pdm_theta = None
        period = d['period']
        condition = (period < 2) & (pdm_theta < 0.8)

        if condition:
            with open(logpath, 'w') as f:
                msg = (
                    f"{starid}: got PDMtheta={pdm_theta:.3f} and "
                    f"P={period:.3f}; exitcode2."
                )
                print(msg)
                f.writelines(msg)
            continue

        # if we have a sufficient periodic signal, then count phased local
        # minima, after normalizing and smoothing.
        cd = {
            'method': 'psplinefilt_findpeaks',
            'height': '2_P2P', # could be "5_MAD", 1e-3, etc.
            'binsize_phase_units': 0.01,
            'width': 2,
            'window_length_phase_units': 0.1,
            'max_splines': 10,
            'height_limit': 1e-3,
            'pre_normalize': True
        }
        r = count_phased_local_minima(
            d['times'], d['fluxs'], d['t0'], d['period'],
            method=cd['method'],
            binsize_phase_units=cd['binsize_phase_units'],
            height=cd['height'], width=cd['width'],
            window_length_phase_units=cd['window_length_phase_units'],
            max_splines=cd['max_splines'],
            height_limit=cd['height_limit'],
            pre_normalize=cd['pre_normalize']
        )

        _pklpath = os.path.join(cachedir, f"{starid}_findpeaks_result.pkl")
        if not os.path.exists(_pklpath):
            with open(_pklpath, 'wb') as f:
                pickle.dump(r, f)
                print(f'Made {_pklpath}')

        #
        # if there are >=3 local minima at whatever confidence, make some plots
        #

        if r['N_peaks'] < 3:
            with open(logpath, 'w') as f:
                msg = f"{starid}: got N={r['N_peaks']} peaks; exitcode3."
                print(msg)
                f.writelines(msg)
            continue

        # make the phased plot
        eval_dict = None
        plotdir = cachedir
        plot_dipcountercheck(r, d, eval_dict, plotdir, starid)

        outpath = join(cachedir, f'{starid}_phase.png')
        titlestr = f"TIC{starid}".replace("_", " ")
        if not os.path.exists(outpath):
            plot_phased_light_curve(
                d['times'], d['fluxs'], d['t0'], d['period'], outpath,
                titlestr=titlestr, binsize_minutes=10, findpeaks_result=r
            )
        else:
            print(f"Found {outpath}")

if __name__ == "__main__":
    find_CPVs()
