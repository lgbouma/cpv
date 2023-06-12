"""
This is a simple pipeline for identifying CPVs.

Contents:
    | find_CPV: given ticid, determine if it's a CPV
    | find_CPVs: thin wrapper
"""
#############
## LOGGING ##
#############
import logging
from complexrotators import log_sub, log_fmt, log_date_fmt

LOCAL_DEBUG = 1
DEBUG = False
if DEBUG:
    level = logging.DEBUG
else:
    level = logging.INFO
LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=level,
    style=log_sub,
    format=log_fmt,
    datefmt=log_date_fmt,
    force=True
)

LOGDEBUG = LOGGER.debug
LOGINFO = LOGGER.info
LOGWARNING = LOGGER.warning
LOGERROR = LOGGER.error
LOGEXCEPTION = LOGGER.exception

#############
## IMPORTS ##
#############
import os, pickle
from os.path import join
from glob import glob
import numpy as np, pandas as pd

from complexrotators.paths import LOCALDIR, SPOCDIR, TARGETSDIR
from complexrotators.getters import (
    _get_lcpaths_given_ticid, _get_local_lcpaths_given_ticid,
    _get_lcpaths_fromlightkurve_given_ticid
)
from complexrotators.lcprocessing import (
    cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
)
from complexrotators.plotting import (
    plot_phased_light_curve, plot_dipcountercheck, plot_cpvvetter
)

from complexrotators import pipeline_utils as pu

def get_ticids(sample_id):

    if sample_id == 'debug':
        ticids = [
        "201789285",
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
        #"245902096",
        #"150068381",
        #"177309964",
        #"118769116",
        #"245868207",
        #"245874053",
        #"59129133"
        ]

    elif sample_id == '30pc_mkdwarf':

        df = pd.read_csv(join(SPOCDIR, "gaia_X_spoc2min_merge.csv"))

        sel = (
            (df["M_G"] > 4)
            &
            (df["bp_rp"] > 1.5)
            &
            (df["TESSMAG"] < 16)
            &
            (df["parallax"] > (100/3))
        )

        sdf = df[sel]
        ticids = np.unique(list(sdf["TICID"].astype(str)))

        N_stars_to_search = len(ticids)
        N_lcs_to_search = len(sdf)

    # e.g.,
    #30to50pc_mkdwarf
    #50to60pc_mkdwarf
    #60to70pc_mkdwarf
    #70to85pc_mkdwarf
    elif 'pc_mkdwarf' in sample_id and 'to' in sample_id:

        lower = int(sample_id.split("to")[0])
        upper = int(sample_id.split("to")[1].split("pc")[0])

        df = pd.read_csv(join(SPOCDIR, "gaia_X_spoc2min_merge.csv"))

        sel = (
            (df["M_G"] > 4)
            &
            (df["bp_rp"] > 1.5)
            &
            (df["TESSMAG"] < 16)
            &
            (df["parallax"] <= 1e3*(1/lower))
            &
            (df["parallax"] > 1e3*(1/upper))
        )

        sdf = df[sel]
        ticids = np.unique(list(sdf["TICID"].astype(str)))

        N_stars_to_search = len(ticids)
        N_lcs_to_search = len(sdf)

    elif sample_id == 'rahul_20230501':
        df = pd.read_csv(
            join(TARGETSDIR, '20230501_RAHUL_FULL_LIST_NO_DUPLICATES.csv')
        )
        ticids = np.unique(df.ticid.astype(str))

        N_stars_to_search = len(ticids)
        N_lcs_to_search = -1

    else:
        raise NotImplementedError

    LOGINFO(42*'-')
    LOGINFO(f"{sample_id}")
    LOGINFO(f"N_stars_to_search = {N_stars_to_search}...")
    LOGINFO(f"N_lcs_to_search = {N_lcs_to_search}...")

    return ticids



def find_CPV(ticid, sample_id, forcepdf=0):
    """
    ticid: e.g. "289840928"
    sample_id: e.g., "30to50pc_mkdwarf" (used for cacheing)
    forcepdf: if true, will require the pdf plot to be made, even if the usual
        exit code criteria were not met.

    exit code definitions:
        exitcode 2: means periodogram_condition was not met
            periodogram_condition = (period < 2) & (pdm_theta < 0.9)

        exitcode 3: means not enough peaks were found --
            r['N_peaks'] < 3:
    """

    cachedir = join(LOCALDIR, "cpv_finding")
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    cachedir = join(cachedir, sample_id)
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    minexitcode = -1
    cand_logpaths = glob(join(cachedir, f"tess*00{ticid}-*runstatus.log"))
    foundexitcodes = []
    if len(cand_logpaths) > 0:
        for cand_logpath in cand_logpaths:
            st = pu.load_status(cand_logpath)
            if 'exitcode' in st:
                exitcode = st['exitcode']['exitcode']
                foundexitcodes.append(int(exitcode))
        if len(foundexitcodes) > 0:
            minexitcode = np.nanmin(foundexitcodes)

    MINIMUM_EXITCODE = 2
    if forcepdf:
        MINIMUM_EXITCODE = 1
    #1 if any kind of exit means do not rerun
    #2 if only a periodogram or not enoigh dip exit means dont rerun
    if minexitcode >= MINIMUM_EXITCODE:
        LOGINFO(f"TIC{ticid}: found log for {ticid} with exitcode {minexitcode}. skip.")
        return 1

    #
    # get the light curves for all desired sectors and cadences
    #
    if LOCAL_DEBUG:
        lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid)
    else:
        lcpaths = _get_local_lcpaths_given_ticid(ticid)

    #
    # for each light curve (sector / cadence specific), detrend, ((remove
    # flares)), get the best period, and then phase-fold.
    #
    for lcpath in lcpaths:

        # instantiate the log
        lcpbase = os.path.basename(lcpath).replace(".fits", "")
        logpath = join(cachedir, f'{lcpbase}_runstatus.log')
        if not os.path.exists(logpath):
            lcpd = {
                'lcpath': lcpath,
                'ticid': ticid
            }
            pu.save_status(logpath, 'lcpath', lcpd)
            LOGINFO(f"Made {logpath}")

        st = pu.load_status(logpath)
        if 'exitcode' in st:
            exitcode = st['exitcode']['exitcode']
            if minexitcode >= MINIMUM_EXITCODE:
                LOGINFO(f"{lcpbase}: found exitcode {exitcode}. skip.")
                if not forcepdf:
                    continue

        # get the relevant light curve data
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        # get period, t0, and periodogram (PDM or LombScargle)
        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        # require peak PDM theta statistic < 0.9
        # and peak period < 2 days
        pdm_theta = d['lsp']['bestlspval']
        period = d['period']
        periodogram_condition = (period < 2) & (pdm_theta < 0.9)

        psr = {
            'starid': starid,
            'sector': sector,
            'cadence_sec': cadence_sec,
            'period': period, # fine_lsp period
            'pdm_theta': pdm_theta,
            'periodogram_condition': periodogram_condition,
            't0': d['t0'],
            'nbestperiods': d['lsp']['nbestperiods'],
            'nbestlspvals': d['lsp']['nbestlspvals']
        }
        pu.save_status(logpath, 'cpv_periodsearch_results', psr)
        LOGINFO(f"Updated {logpath} with cpv_periodsearch_results")

        if not periodogram_condition:
            LOGINFO(f"{starid}: Θ={pdm_theta:.3f}, P={period:.3f}; exit.")
            exitcode = {'exitcode': 2}
            pu.save_status(logpath, 'exitcode', exitcode)
            if not forcepdf:
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
        pu.save_status(logpath, 'count_phased_local_minima_options', cd)
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

        # save the non-array output to the log file
        save_keys = ['N_peaks', 'peaks_phaseunits', 'peaks', 'properties',
                     't0', 'period', 'binsize_phase_units', 'height', 'width',
                     'p2p_raw', 'p2p_est', 'a_95_5', 'mad', 'mult_fac',
                     '_y_mean', 'nsplines_total', 'nsplines_singlephase']
        outd = {k:v for k,v in r.items() if k in save_keys}
        for k,v in outd.items():
            if isinstance(v, np.ndarray):
                outd[k] = list(v)
            if isinstance(v, dict):
                for _k, _v in v.items():
                    if isinstance(v, np.ndarray):
                        outd[k][_k] = list(_v)
        pu.save_status(logpath, 'count_phased_local_minima_results', outd)

        # save the array output as well, to a pickle file
        _pklpath = os.path.join(cachedir, f"{starid}_findpeaks_result.pkl")
        if not os.path.exists(_pklpath):
            with open(_pklpath, 'wb') as f:
                pickle.dump(r, f)
                LOGINFO(f'Made {_pklpath}')

        #
        # if there are >=3 local minima at whatever confidence, make some plots
        #

        if r['N_peaks'] < 3:
            LOGINFO(f"{starid}: N_peaks={r['N_peaks']}; exit.")
            exitcode = {'exitcode': 3}
            pu.save_status(logpath, 'exitcode', exitcode)
            if not forcepdf:
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
                titlestr=titlestr, binsize_phase=0.005, findpeaks_result=r
            )
        else:
            LOGINFO(f"Found {outpath}")

        outpath = join(cachedir, f'{starid}_cpvvetter.pdf')
        if not os.path.exists(outpath):
            plot_cpvvetter(
                outpath, lcpath, starid, periodsearch_result=d,
                findpeaks_result=r
            )
        else:
            LOGINFO(f"Found {outpath}")

        LOGINFO(f"{starid}: N_peaks={r['N_peaks']}; finished.")
        exitcode = {'exitcode': 1}
        pu.save_status(logpath, 'exitcode', exitcode)


def main():

    #sample_id = 'debug'
    #ticids = get_ticids(sample_id)
    #ticids = ['39970966']
    #for ticid in ticids:
    #    LOGINFO(42*'-')
    #    LOGINFO(f"Beginning {ticid}...")
    #    find_CPV(ticid, sample_id)

    sample_ids = [
        #'debug'
        #'30pc_mkdwarf',
        #'30to50pc_mkdwarf',
        #'50to60pc_mkdwarf',
        #'60to70pc_mkdwarf'
        #'70to85pc_mkdwarf',
        #'85to95pc_mkdwarf',
        #'95to105pc_mkdwarf',
        #'105to115pc_mkdwarf'
        #'115to135pc_mkdwarf'
        #'115to150pc_mkdwarf',
        'rahul_20230501'
    ]

    forcepdf = True # FIXME true only for specific (N<~100 !) samples

    for sample_id in sample_ids:
        ticids = get_ticids(sample_id)
        for ticid in ticids:
            LOGINFO(42*'-')
            LOGINFO(f"Beginning {ticid}...")
            find_CPV(ticid, sample_id, forcepdf=forcepdf)

    LOGINFO("Finished 🎉🎉🎉")


if __name__ == "__main__":
    main()
