import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
from glob import glob
from os.path import join
from complexrotators.paths import LOCALDIR, RESULTSDIR
from complexrotators.getters import _get_lcpaths_fromlightkurve_given_ticid
from complexrotators.lcprocessing import prepare_cpv_light_curve
from astrobase.lcmath import time_bin_magseries

from aesthetic.plot import set_style, savefig

from astropy.timeseries import LombScargle

############
# get data #
############

def build_4029_mask(times, fluxs, t0, period):

    cycle_number = np.floor( (times-t0) / period)

    phase = (times - t0)/period - cycle_number

    poff = 1

    # i'm identifying phases that have been wrapped.  ->  any negative
    # phases (that i read off the timegroup) need 1 added

    # about two-thirds is usually "out of bounds"
    indip = (phase > -1/3+poff) | (phase < 1/3)
    outofdip = ~indip

    ok_intervals = [
        [(1780, 1860), (0.15, 0.2)], # for BTJD times 1780-1860, allow φ=0.15-0.20
        [(1780, 1860), (-0.2+poff, -0.15+poff)],
        [(2880, 2940), (0.03, 0.15)],
    ]
    for ok_interval in ok_intervals:

        time_range, phase_range = ok_interval

        _sel = np.zeros(len(times))

        outofdip |= (
            (times > min(time_range)) & (times < max(time_range)) &
            (phase > min(phase_range)) & (phase < max(phase_range))
        )


    return outofdip


def test0():
    """
    manually examine flux vs time for building the tic 4029 mask
    """
    ticid = '402980664'
    sample_id = '2023catalog_LGB_RJ_concat' # used for cacheing
    period = 18.5611/24
    t0 = 1791.12
    cachedir = join(LOCALDIR, "cpv_finding", sample_id)

    lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid)

    # stack over all sectors
    times, fluxs = [], []
    for lcpath in np.sort(lcpaths):
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)
        times.append(x_obs)
        fluxs.append(y_flat)
    times = np.hstack(np.array(times, dtype=object).flatten())
    fluxs = np.hstack(np.array(fluxs, dtype=object).flatten())

    bd = time_bin_magseries(times, fluxs, binsize=900, minbinelems=1)
    btimes, bfluxs = bd['binnedtimes'], bd['binnedmags']

    outofdip = build_4029_mask(btimes, bfluxs, t0, period)

    #
    # make the plot
    #
    set_style("science")
    fig, ax = plt.subplots(figsize=(20,3))
    ax.plot(btimes, bfluxs, c='k', lw=0.5, zorder=1)

    from astrobase.lcmath import find_lc_timegroups
    _, _groups = find_lc_timegroups(btimes[outofdip], mingap=0.5/24)
    for _g in _groups:
        ax.plot(btimes[outofdip][_g], bfluxs[outofdip][_g], c='red', lw=0.5, zorder=2)

    cyclemin, cyclemax = 0, 65
    ax.set_xlim([t0+cyclemin*period, t0+cyclemax*period])
    ax.set_ylim([0.95, 1.05])

    outdir = join(RESULTSDIR, "4029_mask")
    if not os.path.exists(outdir): os.mkdir(outdir)

    outpath = join(outdir, f"mask_cycle{str(cyclemin).zfill(4)}-{str(cyclemax).zfill(4)}.png")
    savefig(fig, outpath, dpi=400)

    cyclemin, cyclemax = 247, 316
    ax.set_xlim([t0+cyclemin*period, t0+cyclemax*period])
    outpath = join(outdir, f"mask_cycle{str(cyclemin).zfill(4)}-{str(cyclemax).zfill(4)}.png")
    savefig(fig, outpath, dpi=400)

    cyclemin, cyclemax = 1232, 1265
    ax.set_xlim([t0+cyclemin*period, t0+cyclemax*period])
    outpath = join(outdir, f"mask_cycle{str(cyclemin).zfill(4)}-{str(cyclemax).zfill(4)}.png")
    savefig(fig, outpath, dpi=400)

    cyclemin, cyclemax = 1410, 1482
    ax.set_xlim([t0+cyclemin*period, t0+cyclemax*period])
    outpath = join(outdir, f"mask_cycle{str(cyclemin).zfill(4)}-{str(cyclemax).zfill(4)}.png")
    savefig(fig, outpath, dpi=400)


def test1():
    """
    try fitting lomb scargle periodogram sinusoid model to out of dip data
    """
    ticid = '402980664'
    sample_id = '2023catalog_LGB_RJ_concat' # used for cacheing
    period = 18.5611/24
    t0 = 1791.12
    cachedir = join(LOCALDIR, "cpv_finding", sample_id)

    lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid)

    # stack over all sectors
    times, fluxs = [], []
    for lcpath in np.sort(lcpaths):
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)
        times.append(x_obs)
        fluxs.append(y_flat)
    times = np.hstack(np.array(times, dtype=object).flatten())
    fluxs = np.hstack(np.array(fluxs, dtype=object).flatten())

    bd = time_bin_magseries(times, fluxs, binsize=900, minbinelems=1)
    btimes, bfluxs = bd['binnedtimes'], bd['binnedmags']

    outofdip = build_4029_mask(btimes, bfluxs, t0, period)

    err = 1e-4 * np.ones(len(bfluxs))

    # fit lomb scargle model

    nterms = 2
    ls = LombScargle(btimes[outofdip], bfluxs[outofdip], err[outofdip], nterms=nterms)

    period_max = 18.6 / 24
    period_min = 18.5 / 24
    freq, power = ls.autopower(minimum_frequency=1/period_max,
                               maximum_frequency=1/period_min)

    ls_freq_1 = freq[np.argmax(power)]
    ls_period_1 = 1/ls_freq_1

    # note: no need to do on binned times vs non-binned; either is good
    flux_fit_1 = ls.model(btimes, ls_freq_1)
    flux_r = bfluxs - flux_fit_1

    #
    # make the plot
    #
    set_style("science")
    fig, ax = plt.subplots(figsize=(20,3))

    ax.scatter(btimes, bfluxs, c='k', s=0.5, zorder=1)
    ax.plot(btimes, flux_fit_1, c='C0', zorder=2, lw=0.5)

    yoffset = 0.9
    ax.scatter(btimes, flux_r+yoffset, c='k', s=0.5, zorder=1)

    outdir = join(RESULTSDIR, "4029_mask")
    if not os.path.exists(outdir): os.mkdir(outdir)

    cyclemin, cyclemax = 0, 65
    ax.set_xlim([t0+cyclemin*period, t0+cyclemax*period])
    ax.set_ylim([0.8, 1.05])
    outpath = join(outdir, f"lsresid_mask_cycle{str(cyclemin).zfill(4)}-{str(cyclemax).zfill(4)}.png")
    savefig(fig, outpath, dpi=400)

    cyclemin, cyclemax = 247, 316
    ax.set_xlim([t0+cyclemin*period, t0+cyclemax*period])
    outpath = join(outdir, f"lsresid_mask_cycle{str(cyclemin).zfill(4)}-{str(cyclemax).zfill(4)}.png")
    savefig(fig, outpath, dpi=400)

    cyclemin, cyclemax = 1232, 1265
    ax.set_xlim([t0+cyclemin*period, t0+cyclemax*period])
    outpath = join(outdir, f"lsresid_mask_cycle{str(cyclemin).zfill(4)}-{str(cyclemax).zfill(4)}.png")
    savefig(fig, outpath, dpi=400)

    cyclemin, cyclemax = 1410, 1482
    ax.set_xlim([t0+cyclemin*period, t0+cyclemax*period])
    outpath = join(outdir, f"lsresid_mask_cycle{str(cyclemin).zfill(4)}-{str(cyclemax).zfill(4)}.png")
    savefig(fig, outpath, dpi=400)

def get_4029_manual_mask(times, fluxs, cadencenos, get_full=0):

    # made using glue, with the manual_masking_20230617.glu session
    maskpaths = np.sort(
        glob(join(RESULTSDIR, 'tic4029_segments', '*manual_mask*.csv'))
    )

    # out of dip times, fluxs, cadenceno's
    ood_df = pd.concat((pd.read_csv(m) for m in maskpaths))
    ood_df = ood_df.sort_values(by='cadenceno')

    # full out of dip mask
    full_ood = np.in1d(cadencenos, np.array(ood_df.cadenceno))

    # sectors 18 and 19 out of dip
    s18s19_ood = full_ood & (times > 1750) & (times < 1900)

    # sectors 25 and 26 out of dip
    s25s26_ood = full_ood & (times > 1950) & (times < 2100)

    # sectors 53, 58, 59 out of dip
    s53s58s59_ood =  full_ood & (times > 2600)

    if get_full:
        return full_ood

    if get_split:
        return s18s19_ood, s25s26_ood, s53s58s59_ood, full_ood


def test2(nterms=1):
    """
    try fitting lomb scargle periodogram sinusoid model to out of dip data,
    using manual masks and splitting into chunks
    """
    ticid = '402980664'
    sample_id = '2023catalog_LGB_RJ_concat' # used for cacheing
    period = 18.5611/24
    t0 = 1791.12
    cachedir = join(LOCALDIR, "cpv_finding", sample_id)

    lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid)

    # stack over all sectors
    times, fluxs, cadencenos = [], [], []
    for lcpath in np.sort(lcpaths):
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid, cadenceno) = prepare_cpv_light_curve(
             lcpath, cachedir, returncadenceno=1
         )
        times.append(x_obs)
        fluxs.append(y_flat)
        cadencenos.append(cadenceno)
    times = np.hstack(np.array(times, dtype=object).flatten())
    fluxs = np.hstack(np.array(fluxs, dtype=object).flatten())
    cadencenos = np.hstack(np.array(cadencenos, dtype=object).flatten())

    model_id = f"manual_20230617_mask_v0_nterms{nterms}"
    full_ood = get_4029_manual_mask(times, fluxs, cadencenos, get_full=1)

    # btjd
    timechunks = [
        (1750, 1900), # s18/s19
        (1950, 2100), # s25/s26
        (2740, 2770), # s53
        (2880, 2940), # s58/s59
    ]

    outtimes, outfluxs = [], []
    for ix, timechunk in enumerate(timechunks):

        # this time chunk
        sel = (times > timechunk[0]) & (times < timechunk[1])
        # out of dip during this time chunk
        sel_ood = full_ood & (times > timechunk[0]) & (times < timechunk[1])

        # fit lomb scargle model
        ls = LombScargle(times[sel_ood], fluxs[sel_ood],
                         1e-4*np.ones_like(fluxs[sel_ood]), nterms=nterms)

        period_max = 18.6 / 24
        period_min = 18.5 / 24
        freq, power = ls.autopower(minimum_frequency=1/period_max,
                                   maximum_frequency=1/period_min)

        ls_freq_1 = freq[np.argmax(power)]
        ls_period_1 = 1/ls_freq_1

        # note: no need to do on binned times vs non-binned; either is good
        flux_fit_1 = ls.model(times[sel], ls_freq_1)
        flux_r = fluxs[sel] - flux_fit_1

        #
        # make the plot
        #
        plt.close("all")
        set_style("science")
        fig, ax = plt.subplots(figsize=(25,3))

        ax.scatter(times[sel], fluxs[sel], c='k', s=0.5, zorder=1, linewidths=0)
        ax.plot(times[sel], flux_fit_1, c='C0', zorder=2, lw=0.5)

        yoffset = 0.9
        ax.scatter(times[sel], flux_r+yoffset, c='k', s=0.5, zorder=1, linewidths=0)

        outtimes.append(times[sel])
        outfluxs.append(flux_r)

        outdir = join(RESULTSDIR, "4029_mask")
        if not os.path.exists(outdir): os.mkdir(outdir)

        if ix == 0:
            cyclemin, cyclemax = 0, 65
        elif ix == 1:
            cyclemin, cyclemax = 247, 316
        elif ix == 2:
            cyclemin, cyclemax = 1232, 1265
        elif ix == 3:
            cyclemin, cyclemax = 1410, 1482

        ax.set_xlim([t0+cyclemin*period, t0+cyclemax*period])
        ax.set_ylim([0.8, 1.05])
        outpath = join(outdir, f"lsresid_{model_id}_cycle{str(cyclemin).zfill(4)}-{str(cyclemax).zfill(4)}.png")
        savefig(fig, outpath, dpi=400)

    outpath = join(outdir, f"lc_lsresid_{model_id}.csv")

    times = np.hstack(np.array(outtimes, dtype=object).flatten())
    r_fluxs = np.hstack(np.array(outfluxs, dtype=object).flatten())

    out_df = pd.DataFrame({
        'time': times,
        'r_flux': r_fluxs # resdu
    })

    out_df.to_csv(outpath, index=False)
    print(f"Wrote {outpath}")

if __name__ == "__main__":

    do_test0 = 0
    do_test1 = 0
    do_test2 = 1

    if do_test0:
        test0()
    if do_test1:
        test1()
    if do_test2:
        test2(nterms=1)
        test2(nterms=2)

    # TODO:
    #
    # one idea for smarter automated mask construction...
    # use the derivatives!
    # if |df/dt| is less than some value over a large enough number of points,
    # then you are in a "smooth region".
    # https://derivative.readthedocs.io/en/latest/notebooks/Examples.html#Examples
    #
    # alternatively:
    # * try extra terms in the LS fit
    # * manually fine-tune the mask more
    # * let the LS fit change over the three main chunks (s18/19,  s25/26, s>50)

    pass
