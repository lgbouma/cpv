import os
from astropy.io import fits
from glob import glob
from os.path import join
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from complexrotators.paths import RESULTSDIR
from astropy.time import Time
from astropy import units as u
from numpy import array as nparr
from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator
)


from cdips.utils.lcutils import astropy_utc_time_to_bjd_tdb
from aesthetic.plot import set_style, savefig

def plot_ew_timeseries(
    linename = 'HGamma',
    linekey = 'Hγ',
    utcdatestr = '20231112',
    inst = 'DBSP',
    ra = 16.733175*u.deg,
    dec = 80.45945*u.deg
):

    # parse input
    if inst == 'HIRES' and linekey == 'Hβ':
        return None, None, None
    if inst == 'HIRES' and linekey in ['Hα', 'Hγ', 'Hδ']:
        if utcdatestr == '20231123':
            hiresdate = 'j531'
        else:
            raise NotImplementedError(
                'manually input utc utcdatestr-> hires date mapping'
            )
        chip = 'b'
        if linekey == 'Hα':
            chip = 'i'

    dbsp_bluelines = ['Hβ', 'Hγ']
    dbsp_redlines = ['Hα']
    if inst == 'DBSP' and linekey in dbsp_bluelines:
        side = 'blue'
    elif inst == 'DBSP' and linekey in dbsp_redlines:
        side = 'red'

    if inst == 'DBSP':
        fitsdir = f'/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX/{utcdatestr}/p200_dbsp_{side}_A/Science'
        ewdir = f'/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX/{utcdatestr}/p200_dbsp_{side}_A/Balmer_EWs'
        fitspaths = np.sort(glob(join(fitsdir, f"spec1d_{side}*LP_12-502*fits")))
    elif inst == 'HIRES':
        fitsdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC402980664_RDX'
        ewdir = '/Users/luke/Dropbox/proj/cpv/results/HIRES_results/Balmer_EWs'
        fitspaths = np.sort(glob(join(fitsdir, f"{chip}{hiresdate}*fits")))

    print(len(fitspaths))

    times, keys = [], []
    for ix, f in enumerate(fitspaths):
        hl = fits.open(f)

        # DBSP: shutter open time, from "UTSHUT" keyword
        # HIRES: also shutter open time
        mjd = np.float64(hl[0].header['MJD'])
        # exposure time in days
        texp = (hl[0].header['EXPTIME']*u.second).to(u.day).value
        # exposure midtime (mjd)
        time = mjd + texp
        key = os.path.basename(f).replace(".fits","")
        times.append(time)
        keys.append(key)

        hl.close()

    csvpaths = [glob(join(ewdir, '*'+linename+"*"+key+'*_results.csv'))[0] for key in keys]

    dfs = [pd.read_csv(f) for f in csvpaths]

    fitted_ews = [df['Fitted_EW_mA'].iloc[0] for df in dfs] # milliangstr
    perr = [df['Fitted_EW_mA_perr'].iloc[0] for df in dfs]
    merr = [df['Fitted_EW_mA_merr'].iloc[0] for df in dfs]
    errs = np.array([merr,perr])

    times = np.array(times)+2400000.5 # mjd to jd

    t = Time(times, format='jd', scale='utc')

    t_bjd_tdb = astropy_utc_time_to_bjd_tdb(
        t, ra, dec, observatory='earthcenter', get_barycorr=0
    )

    set_style('clean')
    fig, ax = plt.subplots()
    ax.errorbar(t_bjd_tdb - 2460255, fitted_ews, yerr=0.75*errs)
    ax.set_ylabel(f'{linekey} EW [mA]')
    ax.set_xlabel('BJD_TDB - 2460255')
    outdir = join(RESULTSDIR, 'EW_results')
    if not os.path.exists(outdir): os.mkdir(outdir)
    outpath = join(outdir, f'{inst}_{linename}_vs_time_{utcdatestr}.png')
    savefig(fig, outpath)

    ptimes, pews, perrs = t_bjd_tdb - 2460255, fitted_ews, 0.75*errs
    return ptimes, pews, perrs


def t_to_phase(t, dofloor=0):
    t0 = 2450000 + 1791.12
    period = 18.5611/24
    φ = (t-t0)/period
    if dofloor:
        φ = (t-t0)/period- np.floor( (t-t0)/period )
    return φ

def phase_to_t(φ):
    # note: kind of wrong... there is no inverse...
    t0 = 2450000 + 1791.12
    period = 18.5611/24
    return φ * period + t0

def plot_stack_ew_timeseries(ha, hb, hc, utcdatestr, uinsts):

    linekeys = 'Hα,Hβ,Hγ'.split(',')
    ewinfo = [ha, hb, hc]

    set_style('clean')
    fig, axs = plt.subplots(figsize=(3,6), nrows=len(ewinfo), sharex=True)

    for ix, (ax, ewi, lk) in enumerate(zip(axs, ewinfo, linekeys)):

        ptimes, pews, perrs = ewi

        ax.errorbar(ptimes, nparr(pews)/(1e3), yerr=nparr(perrs)/(1e3), c='k',
                    ecolor='k', elinewidth=0.5)

        ax.set_ylabel(f'{lk} EW [$\AA$]')

        ax2 = ax.twiny()
        t = ptimes + 2460255 # convert to bjdtdb
        phases = t_to_phase(t) - np.min(t_to_phase(t)).astype(int)
        if max(phases) > 1:
            phases -= 1
        ax2.errorbar(phases, nparr(pews)/(1e3),
                     yerr=nparr(perrs)/(1e3), alpha=0)
        if ix == 0:
            ax2.set_xlabel("Phase, φ")

    axs[-1].set_xlabel('BJD_TDB - 2460255')

    fig.tight_layout()

    outdir = join(RESULTSDIR, 'EW_results')
    if not os.path.exists(outdir): os.mkdir(outdir)
    outpath = join(outdir, f'{uinsts}_stackew_vs_time_DBSP_UT{utcdatestr}.png')
    savefig(fig, outpath)

def plot_stack_ew_vs_phase(has, hbs, hcs, utcdatestrs, insts):

    linekeys = 'f,Hα,Hβ,Hγ'.split(',')

    datecolors = [f'C{ix}' for ix in range(len(utcdatestrs))]

    set_style('clean')
    fig, axs = plt.subplots(figsize=(3,7), nrows=len(linekeys), sharex=True)

    # iterate over dates / insts
    for _id, (ha,hb,hc,utcdatestr,datec,i) in enumerate(
        zip(has, hbs, hcs, utcdatestrs, datecolors, insts)
    ):

        ewinfo = [None, ha, hb, hc]

        axs[0].text(0.97,0.03+_id*0.06, utcdatestr,
                    transform=axs[0].transAxes, ha='right',va='bottom',
                    color=datec)

        # iterate over lines
        for ix, (ax, ewi, lk) in enumerate(zip(axs, ewinfo, linekeys)):

            if i == 'HIRES' and lk == 'Hβ':
                continue

            if lk != 'f':
                ptimes, pews, perrs = ewi
                t = ptimes + 2460255 # convert to bjdtdb
                phase = t_to_phase(t, dofloor=1)

                ax.scatter(phase, nparr(pews)/(1e3), c=datec, alpha=0.9,
                           linewidths=0, zorder=2, s=3)
                ax.errorbar(phase, nparr(pews)/(1e3), yerr=nparr(perrs)/(1e3),
                            c=datec, ecolor=datec, elinewidth=0.5, lw=0, alpha=0.5,
                            zorder=1)

                ax.set_ylabel(f'{lk} EW [$\AA$]')
                if lk == 'Hα':
                    ax.yaxis.set_major_locator(MultipleLocator(1))
                else:
                    ax.yaxis.set_major_locator(MultipleLocator(2))

            else:
                TIERRASDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/TIERRAS'
                if utcdatestr == '20231111':
                    df = pd.read_csv(
                        join(TIERRASDIR, "20231110_TIC402980664_circular_fixed_ap_phot_13.csv")
                    )
                elif utcdatestr == '20231112':
                    df = pd.read_csv(
                        join(TIERRASDIR, "20231110_TIC402980664_circular_fixed_ap_phot_21.csv")
                    )
                else:
                    continue
                t = np.array(df['BJD TDB'])
                phase = t_to_phase(t, dofloor=1)
                flux = np.array(df['Target Relative Flux'])
                flux_err = np.array(df['Target Relative Flux Error'])

                x0 = 73
                ax.scatter(phase, 1e3*flux - x0, c=datec, alpha=0.9,
                           linewidths=0, zorder=2, s=2)
                ax.errorbar(phase, 1e3*flux - x0, yerr=1e3*flux_err,
                            c=datec, ecolor=datec, elinewidth=0.5, lw=0, alpha=0.5,
                            zorder=1)

                ax.set_ylabel('$f_{\mathrm{TIERRAS}}$ [%]')
                ax.set_ylim([70 - x0,76 - x0])


    axs[-1].set_xlabel("Phase, φ")

    fig.tight_layout(h_pad=0)

    outdir = join(RESULTSDIR, 'EW_results')
    if not os.path.exists(outdir): os.mkdir(outdir)
    uinsts = "_".join(np.unique(insts))
    outpath = join(outdir, f'stackew_vs_phase_{uinsts}_{"_".join(utcdatestrs)}.png')
    savefig(fig, outpath)



if __name__ == "__main__":

    utcdatestrs = "20231111,20231112,20231123".split(",")
    insts = "DBSP,DBSP,HIRES".split(",")
    uinsts = "_".join(np.unique(insts))

    has, hbs, hcs = [], [], []
    for d,i in zip(utcdatestrs, insts):
        ha = plot_ew_timeseries(
            linename='HAlpha', linekey='Hα', utcdatestr=d, inst=i
        )
        if i == 'DBSP':
            hb = plot_ew_timeseries(
                linename='HBeta', linekey='Hβ', utcdatestr=d, inst=i
            )
        hc = plot_ew_timeseries(
            linename='HGamma', linekey='Hγ', utcdatestr=d, inst=i
        )
        plot_stack_ew_timeseries(ha, hb, hc, d, uinsts=uinsts)
        has.append(ha)
        hbs.append(hb)
        hcs.append(hc)

    plot_stack_ew_vs_phase(has, hbs, hcs, utcdatestrs, insts)
