import os
from astropy.io import fits
from glob import glob
from os.path import join
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from complexrotators.paths import RESULTSDIR
from astropy.time import Time
from astropy import units as u, constants as const
from numpy import array as nparr
from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator
)


from cdips.utils.lcutils import astropy_utc_time_to_bjd_tdb
from aesthetic.plot import set_style, savefig

UTCDICT = {
    # LP 12-502
    '20231111': '20231111: DBSP + TIERRAS',
    '20231112': '20231112: DBSP + TIERRAS',
    '20231123': '20231123: HIRES',
    '20231203': '20231203: HIRES + TIERRAS',
    '20231207': '20231207: DBSP + KeplerCam',
    '20240115': '20240115: DBSP + FlareCam',
    '20241104': '20241104: DBSP + TESS',
    '20241105': '20241105: DBSP + TESS',
    '20241106': '20241106: DBSP + TESS',
    '20241107': '20241107: DBSP + TESS',
    '20241108': '20241108: DBSP + TESS',
    # DG CVn
    '20240428': '20240428: DBSP + TIERRAS',
}

def plot_ew_timeseries(
    linename = 'HGamma',
    linekey = 'Hγ',
    utcdatestr = '20231112',
    inst = 'DBSP',
    photinst = 'tess',
    targetid = None,
    ra = None,
    dec = None,
    style = None,
):

    assert isinstance(targetid, str)
    assert not (ra is None)
    assert not (dec is None)
    if style is None:
        style = 'clean'

    # parse input
    if inst == 'HIRES' and linekey == 'Hβ':
        return None, None, None
    if inst == 'HIRES' and linekey in ['Hα', 'Hγ', 'Hδ']:
        if utcdatestr == '20231123':
            hiresdate = 'j531'
        elif utcdatestr == '20231203':
            hiresdate = 'j533'
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
        dbdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX/'
        fitsdir = join(dbdir, f'{utcdatestr}/p200_dbsp_{side}_A/Science')
        ewdir = join(dbdir, f'{utcdatestr}/p200_dbsp_{side}_A/Balmer_EWs')
        fitspaths = np.sort(glob(join(fitsdir, f"spec1d_{side}*{targetid}*fits")))

    elif inst == 'HIRES':
        hrdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/'
        fitsdir = join(hrdir, 'TIC402980664_RDX')
        ewdir = '/Users/luke/Dropbox/proj/cpv/results/HIRES_results/Balmer_EWs'
        fitspaths = np.sort(glob(join(fitsdir, f"{chip}{hiresdate}*fits")))

    print(len(fitspaths))
    assert len(fitspaths) > 0

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
    speccsvpaths = [glob(join(ewdir, '*'+linename+"*"+key+'*_fittedslice.csv'))[0] for key in keys]

    dfs = [pd.read_csv(f) for f in csvpaths]

    # NOTE: The "Fitted EWs" seem janky!
    # fitted_ews = [df['Fitted_EW_mA'].iloc[0] for df in dfs] # milliangstr
    fitted_ews = [np.abs(df['EW_mA'].iloc[0]) for df in dfs] # milliangstr
    perr = [df['Fitted_EW_mA_perr'].iloc[0] for df in dfs]
    merr = [df['Fitted_EW_mA_merr'].iloc[0] for df in dfs]
    errs = np.array([merr,perr])

    times = np.array(times)+2400000.5 # mjd to jd

    t = Time(times, format='jd', scale='utc')

    t_bjd_tdb = astropy_utc_time_to_bjd_tdb(
        t, ra, dec, observatory='earthcenter', get_barycorr=0
    )

    ##########################################
    ##########################################
    set_style(style)
    fig, ax = plt.subplots()
    if targetid == 'LP_12-502':
        t_offset = 2460255 # 0 #min(t_bjd_tdb) # 2460255
    elif targetid == 'DG_CVn':
        t_offset = 0 # int(np.floor(min(t_bjd_tdb)))
    else:
        raise NotImplementedError

    ax.errorbar(t_bjd_tdb - t_offset, fitted_ews, yerr=0.75*errs)
    ax.set_ylabel(f'{linekey} EW [m$\AA$]')
    ax.set_xlabel('BJD$_{\mathrm{TDB}}$ - '+f'{t_offset}')
    #ax.set_xlabel('Time [hours]')

    outdir = join(RESULTSDIR, 'EW_results')
    if not os.path.exists(outdir): os.mkdir(outdir)
    outpath = join(
        outdir,
        f'{targetid}_{inst}_{linename}_vs_time_{utcdatestr}_{style}.png'
    )
    savefig(fig, outpath)

    # bonus: phot overlap
    if utcdatestr in [
        '20240428', # tierras
        '20241104', '20241105', '20241106', '20241107', '20241108' # lp12-502
    ]:

        plt.close("all")

        # get tess data
        if photinst == 'tess':
            TESSDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/tess'
            fitspath = join(
                TESSDIR,
                'tess2024300212641-s0085-0000000402980664-0282-s_lc.fits'
            )
            hdul = fits.open(fitspath)
            d = hdul[1].data
            time = d['TIME'] + 2457000 # btjd -> bjd_tdb
            flux = d['PDCSAP_FLUX']
            flux_err = d['PDCSAP_FLUX_ERR']
            med_flux = np.nanmedian(flux)
            flux /= med_flux
            flux_err /= med_flux
            ms = 1
        elif photinst == 'tierras':
            TIERRASDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/TIERRAS'
            KEPLERCAMDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/KeplerCam'
            TESSDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/tess'
            if utcdatestr == '20231111':
                inst == 'tierras'
                df = pd.read_csv(
                    join(TIERRASDIR, "20231111_TIC402980664_circular_fixed_ap_phot_13.csv")
                )
            elif utcdatestr == '20231112':
                inst == 'tierras'
                df = pd.read_csv(
                    join(TIERRASDIR, "20231112_TIC402980664_circular_fixed_ap_phot_21.csv")
                )
            elif utcdatestr == '20231203':
                inst == 'tierras'
                df = pd.read_csv(
                    join(TIERRASDIR, "20231203_TIC402980664_circular_fixed_ap_phot_18.csv")
                )
            elif utcdatestr == '20231207':
                inst == 'keplercam'
                df = pd.read_csv(
                    join(KEPLERCAMDIR, "LP12-502_20231208_KeplerCam_g.dat"),
                    delim_whitespace=True
                )
                df['BJD TDB'] = df['BJD_TDB_B']
                df['Target Relative Flux'] = df['rel_flux_T1_n']
                df['Target Relative Flux Error'] = df['rel_flux_err_T1_n']
            elif utcdatestr == '20240428':
                df = pd.read_csv(
                    join(TIERRASDIR, "TIC368129164_global_lc.csv"),
                    comment='#'
                )
                sel = (
                    (df['BJD TDB'] > 2460428.81) &
                    (df['BJD TDB'] < 2460429.2)
                )
                df = df[sel]
                df = df.rename({
                    'Flux': 'Target Relative Flux',
                    'Flux Error': 'Target Relative Flux Error'
                }, axis='columns')

            assert len(df) > 0
            time = np.array(df['BJD TDB'])
            flux = np.array(df['Target Relative Flux'])
            flux_err = np.array(df['Target Relative Flux Error'])
            ms = 1

        set_style(style)
        fig, axs = plt.subplots(figsize=(2,2.5), nrows=2, sharex=1)

        ax = axs[0]
        x0 = np.nanmin(t_bjd_tdb - 2457000)
        c = 'k' if 'wob' not in style else 'white'
        ax.errorbar(
            24*(t_bjd_tdb - 2457000 - x0), 1e-3*np.array(fitted_ews),
            yerr=1e-3*0.75*np.array(errs),
            fmt='o', color=c, ecolor=c,
            elinewidth=0.5, capsize=0, capthick=0, markersize=1
        )
        ax.set_ylabel(f'{linekey} EW '+'[$\mathrm{\AA}$]')
        xmin, xmax = ax.get_xlim()

        ax = axs[1]
        #fn2 = lambda y: 100*(y - np.nanmedian(y))
        ax.errorbar(24*(time - 2457000 - x0),
                    100*(flux-np.nanmedian(flux)),
                    yerr=100*flux_err,
                    fmt='o', color=c, ecolor=c, mfc=c,
                    elinewidth=0.2, capsize=0, capthick=0, markersize=ms, mew=0)
        ax.set_ylabel(f'Δ Flux [%]')
        ax.set_xlabel(f'Time [hours]')
        ax.set_xlim((xmin, xmax))
        ax.set_ylim((-3.1, 3.1))

        outdir = join(RESULTSDIR, 'EW_results')
        if not os.path.exists(outdir): os.mkdir(outdir)
        outpath = join(
            outdir,
            f'{targetid}_{inst}_{linename}_vs_time_{utcdatestr}_with{photinst}_{style}.png'
        )
        fig.tight_layout()
        savefig(fig, outpath)

    ptimes, pews, perrs = t_bjd_tdb - t_offset, fitted_ews, 0.75*errs
    return ptimes, pews, perrs, keys, csvpaths, speccsvpaths


def t_to_phase(t, targetid, dofloor=0):

    if targetid == 'LP_12-502':
        t0 = 2450000 + 1791.12 - 0.3*18.5611/24
        period = 18.5611/24
    elif targetid == 'DG_CVn':
        t0 = 2460428
        period = 6.44/24

    φ = (t-t0)/period
    if dofloor:
        φ = (t-t0)/period- np.floor( (t-t0)/period )

    return φ

def phase_to_t(φ, targetid):

    # note: kind of wrong... there is no inverse...
    if targetid == 'LP_12-502':
        t0 = 2450000 + 1791.12
        period = 18.5611/24
    elif targetid == 'DG_CVn':
        t0 = 2460428
        period = 6.44/24

    return φ * period + t0

def plot_stack_ew_timeseries(ha, hb, hc, utcdatestr, uinsts, targetid):

    linekeys = 'Hα,Hβ,Hγ'.split(',')
    ewinfo = [ha, hb, hc]

    set_style('clean')
    fig, axs = plt.subplots(figsize=(3,6), nrows=len(ewinfo), sharex=True)

    for ix, (ax, ewi, lk) in enumerate(zip(axs, ewinfo, linekeys)):

        ptimes, pews, perrs, __keys, __rescsvpaths, __fittedcsvpaths = ewi

        ax.errorbar(ptimes, nparr(pews)/(1e3), yerr=nparr(perrs)/(1e3), c='k',
                    ecolor='k', elinewidth=0.5)

        ax.set_ylabel(f'{lk} EW [$\AA$]')

        ax2 = ax.twiny()

        if targetid == 'LP_12-502':
            t_offset = 2460255
        elif targetid == 'DG_CVn':
            t_offset = 2460428
        else:
            raise NotImplementedError

        t = ptimes + t_offset # convert back to bjdtdb
        phases = t_to_phase(t, targetid, dofloor=1)

        ax2.errorbar(phases, nparr(pews)/(1e3),
                     yerr=nparr(perrs)/(1e3), alpha=0)
        if ix == 0:
            ax2.set_xlabel("Phase, φ")

    axs[-1].set_xlabel(f'BJD_TDB - {t_offset}')

    fig.tight_layout()

    outdir = join(RESULTSDIR, 'EW_results')
    if not os.path.exists(outdir): os.mkdir(outdir)
    outpath = join(outdir, f'{targetid}_{uinsts}_stackew_vs_time_DBSP_UT{utcdatestr}.png')
    savefig(fig, outpath)

def plot_stack_ew_vs_phase(has, hbs, hcs, utcdatestrs, insts, targetid):

    linekeys = 'f,Hα,Hβ,Hγ'.split(',')

    datecolors = [f'C{ix}' for ix in range(len(utcdatestrs))]

    set_style('clean')
    fig, axs = plt.subplots(figsize=(3,7), nrows=len(linekeys), sharex=True)

    # iterate over dates / insts
    for _id, (ha,hb,hc,utcdatestr,datec,i) in enumerate(
        zip(has, hbs, hcs, utcdatestrs, datecolors, insts)
    ):

        ewinfo = [None, ha, hb, hc]

        axs[0].text(0.97,0.03+_id*0.06, UTCDICT[utcdatestr],
                    transform=axs[0].transAxes, ha='right',va='bottom',
                    color=datec, fontsize='x-small')

        # iterate over lines
        for ix, (ax, ewi, lk) in enumerate(zip(axs, ewinfo, linekeys)):

            if i == 'HIRES' and lk == 'Hβ':
                continue

            if lk != 'f':
                ptimes, pews, perrs, __keys, __rescsvpaths, __fittedcsvpaths = ewi

                if targetid == 'LP_12-502':
                    t_offset = 2460255 # set in plot_ew_timeseries
                elif targetid == 'DG_CVn':
                    t_offset = 2460428
                t = ptimes + t_offset # convert to bjdtdb
                phase = t_to_phase(t, targetid, dofloor=1)

                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                perrs[np.abs(perrs) > 1000] = 500 # one dud cosmic ray...
                pews = nparr(pews)
                pews[pews < 0] = np.nan # one Hγ HIRES point negative
                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
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
                KEPLERCAMDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/KeplerCam'
                TESSDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/tess'
                if utcdatestr == '20231111':
                    inst = 'tierras'
                    df = pd.read_csv(
                        join(TIERRASDIR, "20231111_TIC402980664_circular_fixed_ap_phot_13.csv")
                    )
                elif utcdatestr == '20231112':
                    inst = 'tierras'
                    df = pd.read_csv(
                        join(TIERRASDIR, "20231112_TIC402980664_circular_fixed_ap_phot_21.csv")
                    )
                elif utcdatestr == '20231203':
                    inst = 'tierras'
                    df = pd.read_csv(
                        join(TIERRASDIR, "20231203_TIC402980664_circular_fixed_ap_phot_18.csv")
                    )
                elif utcdatestr == '20231207':
                    inst = 'keplercam'
                    df = pd.read_csv(
                        join(KEPLERCAMDIR, "LP12-502_20231208_KeplerCam_g.dat"),
                        delim_whitespace=True
                    )
                    df['BJD TDB'] = df['BJD_TDB_B']
                    df['Target Relative Flux'] = df['rel_flux_T1_n']
                    df['Target Relative Flux Error'] = df['rel_flux_err_T1_n']
                elif utcdatestr == '20240428':
                    inst = 'tierras'
                    df = pd.read_csv(
                        join(TIERRASDIR, "TIC368129164_global_lc.csv"),
                        comment='#'
                    )
                    sel = (
                        (df['BJD TDB'] > 2460428.81) &
                        (df['BJD TDB'] < 2460429.2)
                    )
                    df = df[sel]
                    df = df.rename({
                        'Flux': 'Target Relative Flux',
                        'Flux Error': 'Target Relative Flux Error'
                    }, axis='columns')
                    assert len(df) > 0
                elif utcdatestr in [
                    '20241104', '20241105', '20241106', '20241107', '20241108'
                ]:
                    inst = 'tess'
                    fitspath = join(
                        TESSDIR,
                        'tess2024300212641-s0085-0000000402980664-0282-s_lc.fits'
                    )
                    hdul = fits.open(fitspath)
                    d = hdul[1].data
                    time = d['TIME'] + 2457000 # btjd -> bjd
                    flux = d['PDCSAP_FLUX']
                    flux_err = d['PDCSAP_FLUX_ERR']
                    med_flux = np.nanmedian(flux)
                    flux /= med_flux
                    flux_err /= med_flux
                    sel = (d['TIME'] > 3618) & (d['TIME'] < 3623)
                    df = pd.DataFrame({
                        'BJD TDB': time[sel],
                        'Target Relative Flux': flux[sel],
                        'Target Relative Flux Error': flux_err[sel]
                    })
                else:
                    continue

                t = np.array(df['BJD TDB'])
                phase = t_to_phase(t, targetid, dofloor=1)
                flux = np.array(df['Target Relative Flux'])
                flux_err = np.array(df['Target Relative Flux Error'])

                if utcdatestr not in ['20231207', '20240428'] and inst=='tierras':
                    x0 = 73
                    yval = 1e3*flux - x0
                    yerr = 1e3*flux_err
                else:
                    yval = 1e2*(flux - 1)
                    yerr = 1e2*flux_err

                ax.scatter(phase, yval, c=datec, alpha=0.9,
                           linewidths=0, zorder=2, s=2)
                ax.errorbar(phase, yval, yerr=yerr,
                            c=datec, ecolor=datec, elinewidth=0.5, lw=0, alpha=0.5,
                            zorder=1)

                #if utcdatestr == '20240428':
                #    phase_end = t_to_phase(2460428.96, targetid, dofloor=1)
                #    ax.vlines(phase_end, -5, 5, colors='darkgray', alpha=0.5,
                #              linestyles=':', zorder=-10, linewidths=0.5)
                #    ax.text(phase_end, 3, '(TIERRAS night end)')

                ax.set_ylabel('$f_{\mathrm{Broadband}}$ [%]')
                ax.set_ylim([-5, 5])


    axs[-1].set_xlabel("Phase, φ")

    fig.tight_layout(h_pad=0)

    outdir = join(RESULTSDIR, 'EW_results')
    if not os.path.exists(outdir): os.mkdir(outdir)
    uinsts = "_".join(np.unique(insts))
    outpath = join(outdir, f'{targetid}_stackew_vs_phase_{uinsts}_{"_".join(utcdatestrs)}.png')
    savefig(fig, outpath)


def plot_movie_stack_ew_vs_phase(targetid, has, hbs, hcs, utcdatestrs, insts,
                                 show_vel=0, noline=0):

    linekeys = 'f,Hα,Hβ,Hγ'.split(',')

    datecolors = [f'C{ix}' for ix in range(len(utcdatestrs))]

    set_style('clean')

    fig = plt.figure(figsize=(7,3))
    axd = fig.subplot_mosaic(
        """
        AABBE
        AACCF
        AADDG
        """#,
        #gridspec_kw={
        #    "width_ratios": [1, 1, 1, 1]
        #}
    )
    axkeys = "A,B,C,D".split(",")

    alltimes, allkeys, allfittedpaths, allphases, allresultspaths = (
        {}, {}, {}, {}, {}
    )
    # iterate over dates / insts
    for _id, (ha,hb,hc,utcdatestr,datec,i) in enumerate(
        zip(has, hbs, hcs, utcdatestrs, datecolors, insts)
    ):

        alltimes[f"{utcdatestr}_{i}"] = {}
        allkeys[f"{utcdatestr}_{i}"] = {}
        allfittedpaths[f"{utcdatestr}_{i}"] = {}
        allresultspaths[f"{utcdatestr}_{i}"] = {}
        allphases[f"{utcdatestr}_{i}"] = {}

        ewinfo = [None, ha, hb, hc]

        axd['A'].text(0.97,0.03+_id*0.06, UTCDICT[utcdatestr],
                    transform=axd['A'].transAxes, ha='right',va='bottom',
                    color=datec, fontsize='x-small')

        # iterate over lines
        for ix, (ax, ewi, lk) in enumerate(zip([axd['A'],axd['B'],axd['C'],axd['D']], ewinfo, linekeys)):

            ax.set_xlim([0.18,1.02])

            if i == 'HIRES' and lk == 'Hβ':
                continue

            if lk != 'f':
                ptimes, pews, perrs, _keys, _resultscsvpaths, _fittedcsvpaths, = ewi

                t = ptimes + 2460255 # convert to bjdtdb
                phase = t_to_phase(t, targetid, dofloor=1)

                alltimes[f"{utcdatestr}_{i}"][lk] = ptimes
                allkeys[f"{utcdatestr}_{i}"][lk] = _keys
                allfittedpaths[f"{utcdatestr}_{i}"][lk] = _fittedcsvpaths
                allresultspaths[f"{utcdatestr}_{i}"][lk] = _resultscsvpaths
                allphases[f"{utcdatestr}_{i}"][lk] = phase

                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                perrs[np.abs(perrs) > 1000] = 500 # one dud cosmic ray...
                pews = nparr(pews)
                pews[pews < 0] = np.nan # one Hγ HIRES point negative
                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                ax.scatter(phase, nparr(pews)/(1e3), c=datec, alpha=0.9,
                           linewidths=0, zorder=2, s=3)
                ax.errorbar(phase, nparr(pews)/(1e3), yerr=nparr(perrs)/(1e3),
                            c=datec, ecolor=datec, elinewidth=0.5, lw=0, alpha=0.5,
                            zorder=1)

                WRITEOUT = 1
                if WRITEOUT:
                    csvname = f'{targetid}_{lk}_{utcdatestr}_{ix}.csv'
                    outdf = pd.DataFrame(
                        {'t':nparr(t),
                         'phase':nparr(phase),
                         'ew':nparr(pews)/(1e3),
                         'ew_merr': nparr(perrs)[0,:]/(1e3),
                         'ew_perr': nparr(perrs)[1,:]/(1e3)
                        }
                    )
                    outdf.to_csv(csvname, index=False)
                    print(f"wrote {csvname}")

                ax.set_ylabel(f'{lk} EW [$\AA$]')
                if lk == 'Hα':
                    ax.yaxis.set_major_locator(MultipleLocator(1))
                else:
                    ax.yaxis.set_major_locator(MultipleLocator(2))

                if lk != 'Hγ':
                    ax.set_xticklabels([])

            else:
                TIERRASDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/TIERRAS'
                KEPLERCAMDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/KeplerCam'
                if utcdatestr == '20231111':
                    df = pd.read_csv(
                        join(TIERRASDIR, "20231111_TIC402980664_circular_fixed_ap_phot_13.csv")
                    )
                elif utcdatestr == '20231112':
                    df = pd.read_csv(
                        join(TIERRASDIR, "20231112_TIC402980664_circular_fixed_ap_phot_21.csv")
                    )
                elif utcdatestr == '20231203':
                    df = pd.read_csv(
                        join(TIERRASDIR, "20231203_TIC402980664_circular_fixed_ap_phot_18.csv")
                    )
                elif utcdatestr == '20231207':
                    df = pd.read_csv(
                        join(KEPLERCAMDIR, "LP12-502_20231208_KeplerCam_g.dat"),
                        delim_whitespace=True
                    )
                    df['BJD TDB'] = df['BJD_TDB_B']
                    df['Target Relative Flux'] = df['rel_flux_T1_n']
                    df['Target Relative Flux Error'] = df['rel_flux_err_T1_n']
                else:
                    continue
                t = np.array(df['BJD TDB'])
                phase = t_to_phase(t, targetid, dofloor=1)
                flux = np.array(df['Target Relative Flux'])
                flux_err = np.array(df['Target Relative Flux Error'])

                if utcdatestr != '20231207':
                    x0 = 73
                    yval = 1e3*flux - x0
                    yerr = 1e3*flux_err
                else:
                    yval = 1e2*(flux - 1)
                    yerr = 1e2*flux_err

                ax.scatter(phase, yval, c=datec, alpha=0.9,
                           linewidths=0, zorder=2, s=2)
                ax.errorbar(phase, yval, yerr=yerr,
                            c=datec, ecolor=datec, elinewidth=0.5, lw=0, alpha=0.5,
                            zorder=1)

                WRITEOUT = 1
                if WRITEOUT:
                    csvname = f'{targetid}_{lk}_{utcdatestr}_{ix}.csv'
                    outdf = pd.DataFrame(
                        {'t':nparr(t),
                         'phase':nparr(phase),
                         'flux':nparr(yval),
                         'flux_err': nparr(yerr)
                        }
                    )
                    outdf.to_csv(csvname, index=False)
                    print(f"wrote {csvname}")

                ax.set_ylabel('Relative flux [%]')
                ax.set_ylim([-5, 5])

    axd['A'].set_xlabel("Phase, φ")
    axd['D'].set_xlabel("Phase, φ")

    fig.tight_layout(h_pad=0, w_pad=0.5)

    outdir = join(RESULTSDIR, 'movie_EW_results')
    if not os.path.exists(outdir): os.mkdir(outdir)
    uinsts = "_".join(np.unique(insts))
    outpath = join(outdir, f'stackew_vs_phase_{uinsts}_{"_".join(utcdatestrs)}.png')
    savefig(fig, outpath)

    time_indices = []
    for utcdatestr, i in zip(utcdatestrs, insts):
        timed = alltimes[f"{utcdatestr}_{i}"]
        phased = allphases[f"{utcdatestr}_{i}"]
        N_times = len(timed['Hα'])
        # by default, take only the times from the red arm...
        for ix, t, phase in zip(range(N_times), timed['Hα'], phased['Hα']):
            time_index = f"{utcdatestr}_{i}_{str(ix).zfill(4)}_{t:.8f}_{phase:.6f}"
            time_indices.append(time_index)

    # iterate over time indices...
    for _time_index in time_indices:

        # TODO FIXME proceed from here...
        plt.close("all")

        fig = plt.figure(figsize=(7,3.5))
        axd = fig.subplot_mosaic(
            """
            AABBE
            AACCF
            AADDG
            """#,
            #gridspec_kw={
            #    "width_ratios": [1, 1, 1, 1]
            #}
        )
        axkeys = "A,B,C,D".split(",")

        _, i, time_index, this_t, this_phase = _time_index.split("_")
        this_phase = float(this_phase)
        #this_t = np.float64(this_t)
        #this_phase = t_to_phase(this_t, dofloor=1)

        for _id, (ha,hb,hc,utcdatestr,datec,i) in enumerate(
            zip(has, hbs, hcs, utcdatestrs, datecolors, insts)
        ):
            #
            # A: do flux vs phase...
            #
            ax = axd['A']

            axd['A'].text(0.97,0.03+_id*0.05, utcdatestr,
                          transform=axd['A'].transAxes, ha='right',va='bottom',
                          color=datec, fontsize='small')

            TIERRASDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/TIERRAS'
            KEPLERCAMDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/KeplerCam'
            if utcdatestr == '20231111':
                df = pd.read_csv(
                    join(TIERRASDIR, "20231111_TIC402980664_circular_fixed_ap_phot_13.csv")
                )
            elif utcdatestr == '20231112':
                df = pd.read_csv(
                    join(TIERRASDIR, "20231112_TIC402980664_circular_fixed_ap_phot_21.csv")
                )
            elif utcdatestr == '20231203':
                df = pd.read_csv(
                    join(TIERRASDIR, "20231203_TIC402980664_circular_fixed_ap_phot_18.csv")
                )
            elif utcdatestr == '20231207':
                df = pd.read_csv(
                    join(KEPLERCAMDIR, "LP12-502_20231208_KeplerCam_g.dat"),
                    delim_whitespace=True
                )
                df['BJD TDB'] = df['BJD_TDB_B']
                df['Target Relative Flux'] = df['rel_flux_T1_n']
                df['Target Relative Flux Error'] = df['rel_flux_err_T1_n']
            else:
                continue
            t = np.array(df['BJD TDB'])
            phase = t_to_phase(t, targetid, dofloor=1)
            flux = np.array(df['Target Relative Flux'])
            flux_err = np.array(df['Target Relative Flux Error'])

            if utcdatestr != '20231207':
                x0 = 73
                yval = 1e3*flux - x0
                yerr = 1e3*flux_err
            else:
                yval = 1e2*(flux - 1)
                yerr = 1e2*flux_err

            ax.scatter(phase, yval, c=datec, alpha=0.9,
                       linewidths=0, zorder=2, s=2)
            ax.errorbar(phase, yval, yerr=yerr,
                        c=datec, ecolor=datec, elinewidth=0.5, lw=0, alpha=0.5,
                        zorder=1)

            if not noline:
                ax.vlines(this_phase, -5, 5, colors='darkgray', alpha=0.5,
                          linestyles='--', zorder=-10, linewidths=0.5)

            ax.set_ylabel('Broadband flux [%]')
            #ax.set_ylim([-5, 5])
            ax.set_ylim([-3.8, 3.8])

            #
            # B, C, D : line EW vs phase
            #

            # iterate over lines
            ewinfo = [None, ha, hb, hc]
            for ix, (ax, ewi, lk) in enumerate(zip([axd['A'],axd['B'],axd['C'],axd['D']], ewinfo, linekeys)):

                if ix == 0:
                    continue
                ax.set_xlim([0.18,1.02])

                if i == 'HIRES' and lk == 'Hβ':
                    continue

                ptimes, pews, perrs, _keys, _csvpaths, _fittedcsvpaths = ewi

                t = ptimes + 2460255 # convert to bjdtdb
                phase = t_to_phase(t, targetid, dofloor=1)

                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                perrs[np.abs(perrs) > 1000] = 500 # one dud cosmic ray...
                pews = nparr(pews)
                pews[pews < 0] = np.nan # one Hγ HIRES point negative
                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                # FIXME HACK CLEANING
                ax.scatter(phase, nparr(pews)/(1e3), c=datec, alpha=0.9,
                           linewidths=0, zorder=2, s=3)
                ax.errorbar(phase, nparr(pews)/(1e3), yerr=nparr(perrs)/(1e3),
                            c=datec, ecolor=datec, elinewidth=0.5, lw=0, alpha=0.5,
                            zorder=1)

                ax.set_ylabel(f'{lk} '+'EW [$\mathrm{\AA}$]')
                if lk == 'Hα':
                    ax.yaxis.set_major_locator(MultipleLocator(1))
                else:
                    ax.yaxis.set_major_locator(MultipleLocator(2))

                if lk != 'Hγ':
                    ax.set_xticklabels([])

                ymin, ymax = ax.get_ylim()
                if not noline:
                    ax.vlines(this_phase, ymin, ymax, colors='darkgray', alpha=0.5,
                              linestyles='--', zorder=-10, linewidths=0.5)
                ax.set_ylim([ymin, ymax])

        axd['A'].set_xlabel("Phase, φ")
        axd['D'].set_xlabel("Phase, φ")

        #
        # E,F,G: spectra at each time!
        #
        __linekeys = 'Hα,Hβ,Hγ'.split(',')
        keywav = {
            'Hα':6562.8,
            'Hβ':4861.35,
            'Hγ':4340.47
        }
        for ax, lk in zip([axd['E'],axd['F'],axd['G']], __linekeys):

            utckey = "_".join(_time_index.split("_")[:2])
            this_index = int(_time_index.split("_")[2])

            fittedcsvpath = allfittedpaths[utckey][lk][this_index]
            resultscsvpath = allresultspaths[utckey][lk][this_index]

            spec_df = pd.read_csv(fittedcsvpath)
            wav,flx = nparr(spec_df.wav), nparr(spec_df.flx)

            res_df = pd.read_csv(resultscsvpath)
            dx = keywav[lk] - res_df.centroid_A.iloc[0]

            # normalize to fifth entry... (!)
            flx /= flx[5]

            if not show_vel:
                ax.plot(wav + dx, flx, c='k', lw=0.5)

            if show_vel:

                if lk == 'Hα':
                    wav0 = 6563#2.8
                elif lk == 'Hβ':
                    wav0 = 4861#.35
                elif lk == 'Hγ':
                    wav0 = 4340#.47

                deltawvlen = ( (wav + dx) - wav0 )
                delta_v = const.c * (deltawvlen / wav0)
                delta_v_kms = delta_v.to(u.km/u.s)

                ax.plot(delta_v_kms, flx, c='k', lw=0.5)

            ax.update({'ylabel': f'{lk} flux'})

            if not show_vel:
                if lk == 'Hα':
                    wav0 = 6563#2.8
                    ax.set_ylim([0.5,5.6])
                    x0 = 2
                    ax.set_xlim([wav0-8, wav0+8])
                    ax.set_xticks([wav0-5, wav0+5])
                elif lk == 'Hβ':
                    wav0 = 4861#.35
                    ax.set_ylim([0.5,5.6])
                    ax.set_xlim([wav0-8, wav0+8])
                    ax.set_xticks([wav0-5, wav0+5])
                elif lk == 'Hγ':
                    ax.set_ylim([0.5,5.6])
                    wav0 = 4340#.47
                    ax.set_xlim([wav0-8, wav0+8])
                    ax.set_xticks([wav0-5, wav0+5])
            else:
                if lk == 'Hα':
                    ax.set_ylim([0.5,5.6])
                elif lk == 'Hβ':
                    ax.set_ylim([0.5,5.6])
                elif lk == 'Hγ':
                    ax.set_ylim([0.5,5.6])

                ax.set_xlim([-500, 500])
                ax.set_xticks([-250, 0, 250])

            ax.set_yticklabels([1,3,5])
            ax.set_yticks([1,3,5])

        axd['G'].set_xlabel('λ [$\mathrm{\AA}$]')
        if show_vel:
            axd['G'].set_xlabel('v [km$\,$s$^{-1}$]')

        #fig.suptitle('LP 12-502, $P$ = 18.6 hr')
        fig.tight_layout(h_pad=0, w_pad=0.5)

        if not noline:
            outdir = join(RESULTSDIR, 'movie_EW_results')
        else:
            outdir = join(RESULTSDIR, 'movie_EW_results_noline')

        if show_vel:
            outdir = join(RESULTSDIR, 'movie_EW_results_showvel')
        if not os.path.exists(outdir): os.mkdir(outdir)
        uinsts = "_".join(np.unique(insts))

        s = ''
        if show_vel:
            s += "_showvel"

        outpath = join(
            outdir,
            f'{_time_index}_stackew_vs_phase_{uinsts}_{"_".join(utcdatestrs)}{s}.png'
        )
        savefig(fig, outpath, writepdf=0)





if __name__ == "__main__":

    # NOTE: if you add new dates, gotta add to UTCDICT

    # utcdatestrs = "20231111,20231112,20231123,20231203,20231207".split(",")
    # insts = "DBSP,DBSP,HIRES,HIRES,DBSP".split(",")

    # utcdatestrs = "20231111,20231112,20231203,20231207".split(",")
    # insts = "DBSP,DBSP,HIRES,DBSP".split(",")

    # utcdatestrs = "20231111,20231112,20231123".split(",")
    # insts = "DBSP,DBSP,HIRES".split(",")

    #utcdatestrs = "20231111,20231112".split(",")
    #insts = "DBSP,DBSP".split(",")

    # utcdatestrs = ["20240115"]
    # insts = ["DBSP"]

    # utcdatestrs = "20231111,20231112,20231207,20240115".split(",")
    # insts = "DBSP,DBSP,DBSP,DBSP".split(",")

    # utcdatestrs = "20231207,20240115".split(",")
    # insts = "DBSP,DBSP".split(",")

    targetid = 'LP_12-502'
    ra = 16.733175*u.deg
    dec = 80.45945*u.deg

    photinst = 'tess'
    utcdatestrs = ["20241108"]
    utcdatestrs = "20241104,20241105,20241106,20241107,20241108".split(",")

    #photinst = 'tierras'
    #utcdatestrs = "20231111,20231112".split(",")

    #targetid = 'DG_CVn'
    #ra = 202.94307059819*u.deg
    #dec = 29.27621104174*u.deg
    #utcdatestrs = ["20240428"]
    #photinst = 'tierras'

    style = 'clean_wob'
    insts = ["DBSP"]*len(utcdatestrs)

    uinsts = "_".join(np.unique(insts))

    if len(utcdatestrs) > 1:  assert len(insts) == len(utcdatestrs)

    has, hbs, hcs = [], [], []
    for d,i in zip(utcdatestrs, insts):
        ha = plot_ew_timeseries(
            linename='HAlpha', linekey='Hα', utcdatestr=d, inst=i,
            targetid=targetid, ra=ra, dec=dec, style=style, photinst=photinst
        )
        if i == 'DBSP':
            hb = plot_ew_timeseries(
                linename='HBeta', linekey='Hβ', utcdatestr=d, inst=i,
                targetid=targetid, ra=ra, dec=dec, style=style, photinst=photinst
            )
        hc = plot_ew_timeseries(
            linename='HGamma', linekey='Hγ', utcdatestr=d, inst=i,
            targetid=targetid, ra=ra, dec=dec, style=style, photinst=photinst
        )
        plot_stack_ew_timeseries(ha, hb, hc, d, uinsts, targetid)
        has.append(ha)
        hbs.append(hb)
        hcs.append(hc)

    plot_stack_ew_vs_phase(has, hbs, hcs, utcdatestrs, insts, targetid)

    if utcdatestrs == "20231111,20231112".split(","):
        plot_movie_stack_ew_vs_phase(
            targetid, has, hbs, hcs, utcdatestrs, insts, show_vel=0, noline=1
        )
        assert 0
        plot_movie_stack_ew_vs_phase(
            targetid, has, hbs, hcs, utcdatestrs, insts, show_vel=1
        )
        plot_movie_stack_ew_vs_phase(
            targetid, has, hbs, hcs, utcdatestrs, insts
        )
