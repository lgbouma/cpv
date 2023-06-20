"""
Exocomet fit lifted from Sebastian Zieba's beta Pic work:
    https://github.com/sebastian-zieba/betaPic_comet

this is mostly scratchwork
"""
import numpy as np, pandas as pd
from astropy import units as u
from astropy import constants as c
import matplotlib.pyplot as plt
import scipy.optimize as op
import os
from os.path import join
import emcee
import corner

from complexrotators.cometmodel import mu_comet_model
from complexrotators.paths import RESULTSDIR
OUTDIR = join(RESULTSDIR, 'test_exocomet_fit')
if not os.path.exists(OUTDIR): os.mkdir(OUTDIR)

from astrobase.lcmath import phase_magseries_with_errs, phase_magseries

def bins(time, flux, binsize):
    n_bins = len(time) // binsize
    indexes = np.array_split(np.arange(len(time)), n_bins)
    print((indexes[0]).size)
    binned_time = np.array([np.mean(time[a]) for a in indexes])
    binned_flux = np.array([np.mean(flux[a]) for a in indexes])
    binned_flux_err = np.array([np.std(flux[a])/np.sqrt(a.size) for a in indexes])

    return (binned_time, binned_flux, binned_flux_err)

def lnlike(theta, t, f, ferr, ulimb, verbose=False):

    tmid, impact, ce, lamb, period = theta

    model = mu_comet_model(t*u.day, f, ferr, tmid*u.day, impact, ce,
                           lamb/u.radian, period*u.day, Ms, Rs, ulimb)

    chi2 = np.power((f-model)/ferr, 2.)

    if verbose:
        print(chi2)

    return -0.5*(np.sum(chi2))

def plotfit(theta, t, f, ferr, ulimb=None):

    tmid, impact, ce, lamb, period = theta

    model = mu_comet_model(t*u.day, f, ferr, tmid*u.day, impact, ce,
                           lamb/u.radian, period*u.day, Ms, Rs, ulimb)

    plt.close("all")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7), sharex=True)

    ax1.errorbar(t, f, yerr=ferr, c='k',fmt='none')

    ax2.set_xlabel('BTJD')
    ax1.set_ylabel('flx')

    ax1.plot(t, model, c='r')
    ax2.errorbar(t, f-model, yerr=ferr, c='k', fmt='none')
    redchisq = np.power((f-model)/ferr,2.)/f.size

    ax2.axhline(y=0, color='k', linestyle='dotted')

    return fig

def plotfit_phase(theta, t, f, ferr, t0, ulimb=None):

    tmid, impact, ce, lamb, period = theta

    model = mu_comet_model(t*u.day, f, ferr, tmid*u.day, impact, ce,
                           lamb/u.radian, period*u.day, Ms, Rs, ulimb)

    plt.close("all")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7), sharex=True)

    _pd = phase_magseries(t, f, period, t0, wrap=0, sort=1)
    ax1.scatter(_pd['phase'], _pd['mags'], c='k', s=1, linewidths=0, zorder=2)

    _pmd = phase_magseries(t, model, period, t0, wrap=0, sort=1)
    ax1.plot(_pmd['phase'], _pmd['mags'], c='r', zorder=1)

    ax2.scatter(_pd['phase'], _pd['mags']-_pmd['mags'], s=1, linewidths=0)

    ax2.set_xlabel('phase')
    ax1.set_ylabel('flx')

    redchisq = np.power((f-model)/ferr,2.)/f.size

    ax2.axhline(y=0, color='k', linestyle='dotted')

    return fig

def main_betapic():

    datadir = os.getcwd()
    datafile = '/betaPic_mag_res.dat'

    times, mags = np.loadtxt(datadir + datafile, usecols=(0, 1)).T
    tmid = 1486.5 * u.day

    tmin = 1485.0
    tmax = 1488.0

    m_zero = -0.0003 # estimate for no comet magnitude of star

    (times_binned, mags_binned, mags_binned_error) = bins(times, mags, 15)

    tmask = (times_binned>tmin) * (times_binned<tmax)
    t = times_binned[tmask]
    f = np.power(10,(mags_binned[tmask]-m_zero)/-2.5)
    ferr = mags_binned_error[tmask]

    fig, ax = plt.subplots(1, 1, figsize=(8, 10))

    ax.errorbar(t, f, yerr=ferr, c='k',fmt='none')
    ax.set_ylim(0.9970, 1.0005)
    ax.set_xlim(1484., 1490)

    ax.set_xlabel('Time - 2457000 [BTJD days]')
    ax.set_ylabel('inst. mag [mag]')

    plt.savefig('secret_data.png')

    #intial parameter 'guesses'
    par = [1.48629045e+03, 4.73633661e-04, 3.67577405e-03, 8.23439812e+02, 2.29640594e+03]
    labels = ["$t_{mid}$", "$b$", "$c_e$", "$\lambda$","$P$"]

    # get an idea of the density compared to the stellar size

    global Ms, Rs
    Rs = 1.53 * c.R_sun # Wang 2016
    Ms = 1.80 * c.M_sun # Wang 2016

    ulimb = 0.275
    # randomly was 0.79 in Sebastian's example

    t_phase, Inorm, rh, lc = mu_comet_model(
        t*u.day, f, ferr, par[0]*u.day, par[1], par[2], par[3]/u.radian,
        par[4]*u.day, Ms, Rs, ulimb, extra=True
    )

    plt.close("all")
    fig, ax = plt.subplots(1,2, figsize=(9,5))
    ax[0].plot(t_phase, Inorm, label='Inorm')
    ax[0].plot(t_phase, rh, label='rh')
    ax[1].plot(t, lc, label='lc')
    ax[0].legend(loc='best', fontsize='xx-small')
    plt.savefig('secret_fig2.png')

    #ce = 0.0045
    #lamb = 1700.
    #period = 3636.

    nll = lambda *args: -lnlike(*args)

    result = op.minimize(nll, par, method='nelder-mead', args=(t,f,ferr,ulimb),
                         options={'maxiter':10000,'xtol': 1e-8, 'disp': True})

    print(result)

    fig = plotfit(result['x'],t,f,ferr, ulimb=ulimb)
    fig.savefig('secret_fit.png')

    del Ms, Rs

def main_tic4029():

    cycle_range = (48,64)
    φ_range = (0.2,0.4)
    period = 18.5611/24
    t0 = 1791.12
    nterms = 2
    model_id = f"manual_20230617_mask_v0_nterms{nterms}"
    ticid = "TIC_402980664"

    # resid flux from subtracting the sinusoid model; made by build_4029_mask.py
    manual_csvpath = f'/Users/luke/Dropbox/proj/cpv/results/4029_mask/lc_lsresid_{model_id}.csv'
    df = pd.read_csv(manual_csvpath)
    t = np.array(df.time)
    f = 1+np.array(df.r_flux)

    from cdips.utils.lcutils import p2p_rms
    ferr = p2p_rms(f)*np.ones(len(f))

    tmin = t0 + period*cycle_range[0]
    tmax = t0 + period*cycle_range[1]

    φmin = φ_range[0]
    φmax = φ_range[1]

    cycle_number = np.floor( (t-t0) / period)
    φ = (t - t0)/period - cycle_number

    tpmask = (t>tmin) & (t<tmax) & (φ>φmin) & (φ<φmax)

    t = t[tpmask]
    f = f[tpmask]
    ferr = ferr[tpmask]

    _pd = phase_magseries_with_errs(t, f, ferr, period, t0, wrap=0, sort=1)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8, 10))

    axs[0].scatter(t, f, c='k', s=1, linewidths=0)
    axs[0].errorbar(max(t)-2, -0.015, yerr=np.nanmedian(ferr), c='k',fmt='none')
    axs[0].set_xlim(tmin, tmax)
    axs[0].set_xlabel('BTJD')
    axs[0].set_ylabel('flux')

    axs[1].scatter(_pd['phase'], _pd['mags'], c='k', s=1, linewidths=0)
    axs[1].errorbar(φmax-0.01, -0.015, yerr=np.nanmedian(_pd['errs']), c='k',fmt='none')
    axs[1].set_xlabel('phase')
    axs[1].set_ylabel('flux')

    for ax in axs:
        ax.set_ylim([0.97, 1.01])

    plt.savefig(join(OUTDIR, f'{ticid}_data.png'))

    # intial parameter 'guesses'
    # t0: obvious
    # c_e: ~ 2x transit depth, in units of %  (2% depth -> 6.5)
    # λ: scale factor, units of inverse radian

    phase_dur_guess = 0.01
    phase_dur_deg = phase_dur_guess*360*u.deg
    phase_dur_rad = phase_dur_deg.to(u.rad).value
    lambda_guess = 1 / phase_dur_rad

    par = [1830.012, 0, 6.5, lambda_guess, period]
    labels = ["$t_{mid}$", "$b$", "$c_e$", "$\lambda$","$P$"]

    global Ms, Rs
    Rs = 0.5 * c.R_sun # total guess
    Ms = 0.25 * c.M_sun # ditto
    ulimb = 0.75 # guess; Claret2000, linear limb darkening, PHOENIX, Teff3000 logg 4-4.5 ish

    t_phase, Inorm, rh, lc = mu_comet_model(
        t*u.day, f, ferr, par[0]*u.day, par[1], par[2], par[3]/u.radian,
        par[4]*u.day, Ms, Rs, ulimb, extra=True
    )

    plt.close("all")
    fig, ax = plt.subplots(1, 3, figsize=(9,5))
    ax[0].plot(t_phase, Inorm/Inorm.max(), label='Inorm')
    ax[0].plot(t_phase, rh/rh.max(), label='rh')
    ax[1].plot(t, lc, label='lc')
    _pd = phase_magseries(t, lc, period, t0, wrap=0, sort=1)
    ax[2].scatter(_pd['phase'], _pd['mags'], c='k', s=1, linewidths=0)
    ax[0].legend(loc='best', fontsize='xx-small')
    fig.tight_layout()
    plt.savefig(join(OUTDIR, f'{ticid}_fig2.png'))

    #FIXME
    # ["$t_{mid}$", "$b$", "$c_e$", "$\lambda$","$P$"]
    theta = 1830.012, 0, 6.5, lambda_guess, period # by eye params
    theta = [1.83001217e+03, 4.76152918e-04, 5.93408535e+00, 1.41942191e+01, 7.73389108e-01] # neldermead
    fig = plotfit_phase(theta, t, f, ferr, t0, ulimb=ulimb)
    fig.savefig(join(OUTDIR, f'{ticid}_manualfit.png'))

    nll = lambda *args: -lnlike(*args)

    #import IPython; IPython.embed()
    #TODO: speed up
    result = op.minimize(nll, par, method='nelder-mead', args=(t,f,ferr,ulimb),
                         options={'maxiter':10000,'xtol': 1e-8, 'disp': True})

    print(result)

    fig = plotfit_phase(result['x'], t, f, ferr, t0, ulimb=ulimb)
    fig.savefig(join(OUTDIR, f'{ticid}_bestfit.png'))

    import IPython; IPython.embed()

    del Ms, Rs


if __name__ == "__main__":

    do_betapic = 0
    do_lp12 = 1

    if do_betapic:
        main_betapic()

    if do_lp12:
        main_tic4029()
