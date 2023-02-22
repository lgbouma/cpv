"""
This is a simple pipeline for identifying CPVs.

Contents:
    find_CPVs: thin wrapper
    find_CPV
"""
import os, pickle, subprocess
from os.path import join
import numpy as np, pandas as pd
from glob import glob

from astropy.io import fits

from complexrotators.paths import LOCALDIR

from complexrotators.lcprocessing import cpv_periodsearch, count_phased_local_minima

from complexrotators.plotting import plot_phased_light_curve

from wotan import flatten

from copy import deepcopy

def find_CPVs():

    # the TICIDs to search
    ticids = [
    #"201789285",
    #"311092148",
    #"332517282",
    #"405910546",
    #"142173958",
    #"300651846",
    #"408188366",
    #"146539195",
    #"177309964",
    #"425933644",
    #"206544316",
    #"224283342",
    "245902096",
    #"150068381",
    #"177309964",
    #"118769116",
    #"245868207",
    #"245874053",
    #"59129133"
    ]

    # TODO: i think you want to implement this routine for dip-counting as its
    # own testable case, with particular assertion statements respected.

    for ticid in ticids:
        find_CPV(ticid)


def prepare_cpv_light_curve(lcpath):

    hl = fits.open(lcpath)
    hdr = hl[0].header
    d = hl[1].data

    # metadata
    sector = hdr["SECTOR"]

    # light curve data
    time = d['TIME']
    flux = d['PDCSAP_FLUX']
    qual = d['QUALITY']

    # remove non-zero quality flags
    sel = (qual == 0)

    x_obs = time[sel]
    y_obs = flux[sel]

    # normalize around 1
    y_obs /= np.nanmedian(y_obs)

    # NOTE: you could consider removing flares using a time-windowed slider
    # here.  however, for purposes of finding the periods, they are a small
    # enough fraction of the duty cycle that they can probably be ignored.

    # what is the cadence?
    cadence_sec = int(np.round(np.nanmedian(np.diff(x_obs))*24*60*60))

    starid = f'{ticid}_S{str(sector).zfill(4)}_{cadence_sec}sec'

    return time, flux, qual, x_obs, y_obs, cadence_sec, sector, star_id



def find_CPV(ticid):

    #
    # get the light curves for all desired sectors and cadences
    #
    SPOCDIR = "/Users/luke/local/SPOCLC"
    lcpaths = glob(join(SPOCDIR, "lightcurves", f"*{ticid}*.fits"))

    if len(lcpaths) == 0:
        p = subprocess.call([
            "scp", f"luke@wh1:/ar1/TESS/SPOCLC/sector*/*{ticid}*.fits",
            join(SPOCDIR, "lightcurves")
        ])

    outdir = join(LOCALDIR, "cpv_finding")
    if not os.path.exists(outdir): os.mkdir(outdir)

    #
    # for each light curve (sector / cadence specific), detrend, remove flares,
    # get the best period, and then phase-fold.
    #
    for lcpath in lcpaths:

        time, flux, qual, x_obs, y_obs, cadence_sec, sector, star_id = (
            prepare_cpv_light_curve(lcpath)
        )

        #
        # "light" detrending by default. (& cache it)
        #
        pklpath = os.path.join(outdir, f"{starid}_dtr_lightcurve.pkl")
        if os.path.exists(pklpath):
            print(f"Found {pklpath}, loading and continuing.")
            with open(pklpath, 'rb') as f:
                lcd = pickle.load(f)
            y_flat = lcd['y_flat']
            y_trend = lcd['y_trend']
            x_trend = lcd['x_trend']
        else:
            y_flat, y_trend = flatten(x_obs, y_obs, window_length=5.0,
                                      return_trend=True, method='median')
            x_trend = deepcopy(x_obs)
            lcd = {'y_flat':y_flat, 'y_trend':y_trend, 'x_trend':x_trend }
            with open(pklpath, 'wb') as f:
                pickle.dump(lcd, f)
                print(f'Made {pklpath}')

        # get t0, period, lsp
        t0 = 'binmin'
        d = cpv_periodsearch(
            x_obs, y_flat, starid, outdir, t0=t0, periodogram_method='pdm'
        )

        # TODO: first, assess period detection significance here!
        # needs a heuristic... based on the overall distribution, i think
        print("NEED TO COME UP WITH A PERIOD DETECTION SIGNIFICANCE CUT HERE")
        import IPython; IPython.embed()

        r = count_phased_local_minima(
            d['times'], d['fluxs'], d['t0'], d['period'],
            binsize_phase_units=0.005, prominence=1e-3, width=3
        )

        _pklpath = os.path.join(outdir, f"{starid}_findpeaks_result.pkl")
        if not os.path.exists(_pklpath):
            with open(_pklpath, 'wb') as f:
                pickle.dump(r, f)
                print(f'Made {_pklpath}')

        #
        # if there are >=3 local minima at whatever confidence, make some plots
        #

        if r['N_peaks'] < 3:
            print(f"{starid}: got N={r['N_peaks']} peaks, finished.")
            continue

        # make the phased plot
        outpath = os.path.join(
            outdir, f'{ticid}_S{str(sector).zfill(4)}_{cadence_sec}sec_phase.png'
        )

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
