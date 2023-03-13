import os
from glob import glob
import pandas as pd, numpy as np
from os.path import join
from complexrotators.paths import LOCALDIR, RESULTSDIR
from complexrotators.plotting import (
    plot_phased_light_curve, plot_dipcountercheck
)
from complexrotators.getters import _get_lcpaths_given_ticid
from complexrotators.lcprocessing import (
    cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
)

from complexrotators import pipeline_utils as pu

def test_dip_counter(ticid, require_classified=True):

    lcpaths = _get_lcpaths_given_ticid(ticid)
    cachedir = join(LOCALDIR, "cpv_finding")
    plotdir = join(RESULTSDIR, "dipcountercheck")
    if not os.path.exists(cachedir): os.mkdir(cachedir)
    if not os.path.exists(plotdir): os.mkdir(plotdir)

    cdf = pd.read_csv("known_cpv_dip_counts.csv")

    for lcpath in np.sort(lcpaths):

        sector = int(os.path.basename(lcpath).split("-")[1][1:].lstrip("0"))

        if require_classified:
            sel = (
                (cdf.ticid.astype(str) == str(ticid))
                &
                (cdf.sector.astype(int) == sector)
            )
            sdf = cdf[sel]

            if len(sdf) == 0:
                print(f'No classifixn for {os.path.basename(lcpath)}. Skip.')
                continue

        # get the relevant light curve data
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        # get period, t0, and periodogram (PDM or LombScargle)
        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        # v4 (splines are better; so is pre-normalizing; so is a height limit
        # that is p2p rms informed)
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

        eval_dict, make_plot = evaluate_dipcounter(starid, r, cd)

        if make_plot:
            plot_dipcountercheck(r, d, eval_dict, plotdir, starid)

        plot_phase = True
        if plot_phase:
            t0 = d['t0']
            #if str(ticid) == "402980664":
            #    if sector in [18, 19]:
            #        t0 = 1791.1372827806442 + 0.21*0.7732752056726647
            outpath = join(plotdir, f"{starid}_phase.png")
            plot_phased_light_curve(
                d['times'], d['fluxs'], t0, d['period'], outpath,
                titlestr=starid.replace("_", " "), binsize_minutes=10,
                xlim=[-0.6,0.6]
            )



def evaluate_dipcounter(starid, r, cd, save_evaldict=True):

    ticid, sstr, cadence = starid.split("_")
    sector = int(sstr[1:].lstrip("0"))

    df = pd.read_csv("known_cpv_dip_counts.csv")

    sel = (
        (df.ticid.astype(str) == ticid)
        &
        (df.sector == sector)
    )

    sdf = df[sel]

    if len(sdf) == 1:
        make_plot = True

    if not len(sdf) == 1:
        print("nothing to evaluate against; will not plot")
        make_plot = False

    # did you get at least the minimum number of dips?
    found_enough_dips = np.all(r['N_peaks'] >= sdf['ndipsmin'])
    found_toofew_dips = np.all(r['N_peaks'] < sdf['ndipsmin'])
    found_toomany_dips = np.all(r['N_peaks'] > sdf['ndipsmax'])

    # did you get the desired range of dip number?
    found_correct_ndips = np.all(
        (r['N_peaks'] >= sdf['ndipsmin'])
        &
        (r['N_peaks'] <= sdf['ndipsmax'])
    )

    class_str = '' if len(sdf) == 0 else sdf['class'].values[0]

    eval_dict = {
        'ticid': ticid,
        'sstr': sstr,
        'cadence': cadence,
        'class': class_str,
        'found_enough_dips': found_enough_dips,
        'found_toofew_dips': found_toofew_dips,
        'found_toomany_dips': found_toomany_dips,
        'found_correct_ndips': found_correct_ndips,
    }
    for k, v in cd.items():
        eval_dict[k] = v
    if make_plot:
        eval_dict['ndips_min'] = int(sdf['ndipsmin'])
        eval_dict['ndips_max'] = int(sdf['ndipsmax'])
    else:
        eval_dict['ndips_min'] = -1
        eval_dict['ndips_max'] = -1

    eval_dict['N_peaks'] = r['N_peaks']
    #eval_dict['peaks_phaseunits'] = r['peaks_phaseunits']
    #eval_dict['props_peak_heights'] = r['properties']['peak_heights']
    eval_dict['height'] = r['height']
    eval_dict['mad'] = r['mad']
    eval_dict['a_95_5'] = r['a_95_5']

    if save_evaldict:
        outdir = join(RESULTSDIR, "dipcountercheck")
        outpath = join(outdir, f'{starid}_dipcountercheck.txt')
        pu.save_status(outpath, 'eval_dict', eval_dict)
        print(f"Wrote {outpath}")

    return eval_dict, make_plot


def evaluate_run_dipcountercheck(ticids):

    outdir = join(RESULTSDIR, "dipcountercheck")
    outpaths = [
        glob(join(outdir, f'{ticid}*_dipcountercheck.txt'))
        for ticid in ticids
    ]
    outpaths = [item for sublist in outpaths for item in sublist]

    rows = []
    for outpath in outpaths:
        _row = pu.load_status(outpath)
        rows.append(pd.DataFrame(dict(_row['eval_dict']), index=[0]))

    df = pd.concat(rows)
    df['imgpath'] = [os.path.basename(p).replace(".txt",".png") for p in outpaths]

    sel = (df.ndips_min.astype(int) >= 0) & (df.ndips_max.astype(int) >= 0)

    sdf = df[sel]

    cols = ['ticid', 'sstr', 'class', 'found_correct_ndips',
            'found_toomany_dips', 'found_toofew_dips', 'imgpath']

    print(sdf[cols].sort_values(by=['class','ticid','sstr']))
    csvpath = join(RESULTSDIR, "dipcountercheck", "run_summary.csv")
    sdf[cols].sort_values(by=['class','ticid','sstr']).to_csv(
        csvpath, index=False
    )
    print(f"Made {csvpath}")

    # require all rotators to have a "correct" number of dips (which in
    # practice means at most two)
    assert np.all(sdf[sdf['class'] == 'rot'].found_correct_ndips == 'True')
    print("All rotators have <= 2 dips.")

    # require all cpvs to have more than the minimum number of dips
    assert np.all(sdf[sdf['class'] == 'cpv'].found_toofew_dips == 'False')
    print("All CPVs have N dips consistent with labels.")


def test_dip_counter_all_stars():

    ticids = [
        ##########################################
        "402980664"
        # # WORKING CASES
        # # normal rotators
        # "294328887"
        # "150068381",
        # "177309964",
        # "149248196", # AB Dor...lotsss of data...
        # "389423271", # Speedy Mic
        # # fav cpvs
        # "206544316",
        # "425933644",
        # "332517282",
        # "300651846",
        # # cpv/eb
        # "146539195",
        # # cpvs that i care less about
        # "201789285",
        # "311092148",
        # "405910546",
        # "142173958",
        # "408188366",
        ## bonus ebs
        # "245834739",
        # "245868207",
        ##########################################
        # cpvs that are fine but i already have enough
        ##"238597707",
        ##"224283342",
        ##"245902096",
        ##"118769116",
        ##"245874053",
        ##"59129133"
    ]

    run_test = 1
    run_eval = 1

    if run_test:
        for ticid in ticids:
            test_dip_counter(ticid)

    if run_eval:
        evaluate_run_dipcountercheck(ticids)


if __name__ == "__main__":
    test_dip_counter_all_stars()
