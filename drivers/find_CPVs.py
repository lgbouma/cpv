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

LOCAL_DEBUG = 0
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
import time as timemod
from os.path import join
from glob import glob
import numpy as np, pandas as pd

from complexrotators.paths import (
    LOCALDIR, SPOCDIR, QLPDIR, TARGETSDIR, TABLEDIR
)
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

def get_ticids(sample_id, lcpipeline):

    if sample_id == 'debug':
        ticids = [
            #"5714469",
            #"219790149"
            #'243499565' # missed, in Sco-Cen from Stauffer2021
            #"57528302" # great TWA disk
            #"234284556"  # tuchor cpv??
            #"407001106" # 3.5d EB or CPV?
            #"260268310" # nice one
            #"359892714" # UCD CPV from 2406.07154, they just plotted poorly
            #'120355394' # yet another odd B star, HD 176582 from Oleg
            #'125843782' # 2M0437
            '219117956'
        ]

        N_stars_to_search = len(ticids)
        N_lcs_to_search = -1

    elif sample_id == 'dovi':
        df = pd.read_csv(
            join(TARGETSDIR, 'Dovi_20240609_ticids.txt')
        )
        get_ticid = lambda x: x.split("_")[0]
        df['ticid'] = df['ticid_sstr'].apply(get_ticid)
        ticids = np.unique(df['ticid'])

        N_stars_to_search = len(ticids)
        N_lcs_to_search = -1

    elif sample_id == 'jan2026_knownCPVs':
        df = pd.read_csv(
            join(
                TABLEDIR, 'jan2026_compilation',
                'concat_R16_S17_S18_B20_S21_Z19_G22_P23_B24_qlp_0to100pc.csv'
            ), sep="|"
        )
        # 167 CPVs
        df = df.drop_duplicates(subset=['ticid'], keep='first')
        ticids = np.unique(list(df["ticid"].astype(str)))

        N_stars_to_search = len(ticids)
        N_lcs_to_search = -1

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

    elif 'pc' in sample_id and 'to' in sample_id and lcpipeline=='qlp':

        lower = int(sample_id.split("to")[0])
        upper = int(sample_id.split("to")[1].split("pc")[0])

        assert upper <= 500 # parsecs

        df = pd.read_csv(join(QLPDIR, "QLP_s1s55_X_GDR2_parallax_gt_2.csv"))

        sel = (
            (df["parallax"] <= 1e3*(1/lower))
            &
            (df["parallax"] > 1e3*(1/upper))
        )

        sdf = df[sel]
        ticids = np.unique(list(sdf["ticid"].astype(str)))

        N_stars_to_search = len(ticids)
        N_lcs_to_search = len(sdf) # this is incorrect; in 2023, it's ~3x N_stars_to_search

    # e.g., 30to50pc_mkdwarf, 50to60pc_mkdwarf, etc.
    elif 'pc_mkdwarf' in sample_id and 'to' in sample_id and lcpipeline=='spoc2min':

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

    elif sample_id == '2023catalog_LGB_RJ_concat':
        df = pd.read_csv(
            join(
                TABLEDIR, '2023_catalog_table',
                '20230613_LGB_RJ_CPV_TABLE_supplemental_selfnapplied_BACKUP.csv'
            ), sep="|"
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



def find_CPV(ticid, sample_id, forcepdf=0, lcpipeline='spoc2min'):
    """
    Args:

    ticid: e.g. "289840928"

    sample_id: e.g., "30to50pc_mkdwarf" (used for cacheing)

    forcepdf: if true, will require the pdf plot to be made, even if the usual
        exit code criteria were not met.

    lcpipeline: "qlp" or "spoc2min"

    exit code definitions:
        exitcode 2: means periodogram_condition was not met
            periodogram_condition = (period < 2) & (pdm_theta < 0.9)

        exitcode 3: means not enough peaks were found --
            r['N_peaks'] < 3:

        exitcode 4: light curve did not have finite values.

        exitcode 5: faulty LC; flatten failed so dip counter failed.

        exitcode 6: insufficient points in LC; dip counter failed.

        exitcode 7: erroneous binning

        exitcode 8: lcpath not found

        exitcode 9: lcpath found but corrupted
    """

    assert lcpipeline in ["qlp", "spoc2min", "spoc2min_tars"]

    cachedir = join(LOCALDIR, "cpv_finding")
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    cachename = f"{lcpipeline}_{sample_id}"
    cachedir = join(cachedir, cachename)
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    minexitcode = -1
    cand_logpaths = glob(join(cachedir, f"*tess*00{ticid}*runstatus.log"))
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
    if minexitcode >= MINIMUM_EXITCODE and not forcepdf:
        LOGINFO(f"TIC{ticid}: found log for {ticid} with exitcode {minexitcode}. skip.")
        return 1

    #
    # get the light curves for all desired sectors and cadences
    #
    if LOCAL_DEBUG:
        lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid, lcpipeline)
    else:
        lcpaths = _get_local_lcpaths_given_ticid(ticid, lcpipeline)

    if sample_id == 'debug':
        #LOGWARNING("Found debug: taking only a single sector.")
        #lcpaths = [np.sort(lcpaths)[0]]
        lcpaths = np.sort(lcpaths)

    #
    # for each light curve (sector / cadence specific), detrend, ((remove
    # flares)), get the best period, and then phase-fold.
    #
    for lcpath in lcpaths:

        if not os.path.exists(lcpath) and lcpipeline=='qlp':
            LOGWARNING(f'Did not find {lcpath}')
            LOGWARNING(f'Trying to curl it directly...')
            # if the light curve does not exist... we will try to curl it...
            searchstr = "/".join(lcpath.split("/")[4:])
            os.popen(
                f"curl -C - -f --create-dirs --output '{lcpath}' 'https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:HLSP/qlp/{searchstr}'"
            )
            # ...and it'll take some time.
            # this is 100% as janky as it looks.  however, unlike the (many)
            # stackoverflow answers that suggest requests, and json, and
            # whatever... this one will work!
            timemod.sleep(10)

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
        try:
            (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
             sector, starid) = prepare_cpv_light_curve(
                 lcpath, cachedir, lcpipeline=lcpipeline
             )
        except FileNotFoundError:
            LOGWARNING(f"{starid}: Failed to find lcpath {lcpath}.")
            exitcode = {'exitcode': 8}
            pu.save_status(logpath, 'exitcode', exitcode)
            continue
        except TypeError:
            LOGWARNING(f"Got corrupt lcpath {lcpath}.")
            exitcode = {'exitcode': 9}
            pu.save_status(logpath, 'exitcode', exitcode)
            continue

        if y_flat is None:
            LOGWARNING(f"{starid}: Failed to get finite light curve.")
            exitcode = {'exitcode': 4}
            pu.save_status(logpath, 'exitcode', exitcode)
            continue

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
        try:
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
        except np.linalg.LinAlgError as e:
            LOGWARNING(f"{starid}: {e} flatten pspline call faulty; SVD did not converge")
            exitcode = {'exitcode': 5}
            pu.save_status(logpath, 'exitcode', exitcode)
            continue
        except ValueError as e:
            LOGWARNING(f"{starid}: {e} insufficient points in light curve")
            exitcode = {'exitcode': 6}
            pu.save_status(logpath, 'exitcode', exitcode)
            continue

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
            try:
                plot_phased_light_curve(
                    d['times'], d['fluxs'], d['t0'], d['period'], outpath,
                    titlestr=titlestr, binsize_phase=0.005, findpeaks_result=r
                )
            except TypeError as e:
                LOGWARNING(f"{starid}: {e} plots failed with erroneous binning.")
                exitcode = {'exitcode': 7}
                pu.save_status(logpath, 'exitcode', exitcode)
                continue

        else:
            LOGINFO(f"Found {outpath}")

        outpath = join(cachedir, f'{starid}_cpvvetter.pdf')
        if not os.path.exists(outpath):
            plot_cpvvetter(
                outpath, lcpath, starid, periodsearch_result=d,
                findpeaks_result=r, lcpipeline=lcpipeline
            )
        else:
            LOGINFO(f"Found {outpath}")

        LOGINFO(f"{starid}: N_peaks={r['N_peaks']}; finished.")
        exitcode = {'exitcode': 1}
        pu.save_status(logpath, 'exitcode', exitcode)


def main():

    #################
    # begin options #
    #################
    forcepdf = 1 # if yes, perhaps also have "LOCALDEBUG" set true..

    # lcpipeline: "qlp", "spoc2min", or "spoc2min_tars"
    lcpipeline = 'spoc2min_tars'

    sample_ids = [
        'jan2026_knownCPVs'
        #'debug',
        #'dovi'
        # ### samples for the next CPV project:
        #'1to20pc'
        # '20to40pc'
        # '40to60pc'
        # '60to80pc'
        # '80to100pc'
        # ### samples for CPV paper #1:
        #'2023catalog_LGB_RJ_concat'
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
        #'rahul_20230501'
    ]

    ###############
    # end options #
    ###############

    for sample_id in sample_ids:

        ticids = get_ticids(sample_id, lcpipeline)

        if len(ticids) > 1000 and forcepdf:
            raise NotImplementedError

        for ticid in ticids:
            LOGINFO(42*'-')
            LOGINFO(f"Beginning {ticid}...")
            find_CPV(ticid, sample_id, forcepdf=forcepdf, lcpipeline=lcpipeline)

    LOGINFO("Finished 🎉🎉🎉")


if __name__ == "__main__":
    main()
