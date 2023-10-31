import numpy as np, matplotlib.pyplot as plt
from os.path import join
import batman

from complexrotators.lcprocessing import (
    cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
)
from complexrotators.plotting import (
    plot_phased_light_curve, plot_dipcountercheck, plot_cpvvetter
)


P0 = 0.5 # give it a 0.5 day fundamental period

def get_ytra_ysec(t, P_eps=0):

    params = batman.TransitParams()
    params.t0 = 0.                       #time of inferior conjunction
    params.per = P0+P_eps                      #orbital period
    params.rp = 0.2                      #planet radius (in units of stellar radii)
    params.a = 3                       #semi-major axis (in units of stellar
    params.inc = 87.                     #orbital inclination (in degrees)
    params.ecc = 0.2                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    params.u = [0.1, 0.3]                #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       #limb darkening model
    params.fp = 0.02
    params.t_secondary = 0.
    m = batman.TransitModel(params, t)    #initializes model
    msec = batman.TransitModel(params, t, transittype="secondary")    #initializes model
    ytra = m.light_curve(params) - 1
    ysec = msec.light_curve(params) - 1

    return ytra, ysec


def get_y0_y1(t):

    P1 = P0*(3/2 + 1e-3)

    phi0 = 2
    A0 = 0.1
    phi1 = 1.12
    A1 = 0.0  # only one sinusoid signal in this example

    y0 = A0 * np.sin(2*np.pi*t/P0 + phi0)
    y1 = A1 * np.sin(2*np.pi*t/P1 + phi1)

    return y0, y1


def main(eps):

    #
    # generate one example
    #
    t = np.linspace(0,28, 20000) # underlying model
    tobs = np.arange(0,28, 2/(60*24)) # tess 120-second cadence

    # P0 = 0.5 days...
    P_eps = eps * P0

    ytra, ysec = get_ytra_ysec(t, P_eps=P_eps)

    y0, _ = get_y0_y1(t)

    ytot = 1 + y0 + ytra

    ytra_obs, ysec_obs = get_ytra_ysec(tobs, P_eps=P_eps)
    y0_obs, _ = get_y0_y1(tobs)

    np.random.seed(42)
    yerr = np.random.normal(loc=0, scale=0.02, size=len(tobs))

    yobs = 1 + y0_obs + ytra_obs + yerr

    #
    # plot the truth, and the observed fake data
    #
    fig, axs = plt.subplots(nrows=3, figsize=(10,6))

    axs[0].plot(t, ytot)
    axs[1].scatter(tobs, yobs, c='k', s=1, linewidths=0)
    ys = [y0+ytra, y0, ytra]
    labels = ['sinusoid plus eclipse', 'sinusoid at P0', 'eclipses at P0+P_eps']
    for ix, (y, l) in enumerate(zip(ys, labels)):
        axs[2].scatter(t, y, c=f'C{ix}', s=1, linewidths=0, label=l)
    axs[2].legend(loc='best', fontsize='small')
    axs[2].set_xlabel('time [days]')

    starid = f"dephasing_eps{eps:.1e}"
    fig.savefig(f'test_output/test_{starid}.png', dpi=400, bbox_inches='tight')
    plt.close("all")

    #
    # process fake observed data to see if we would be fooled
    #
    cachedir = "./test_output/"

    time, flat_flux = tobs*1., yobs*1.

    # NOTE to Saul et al: everything down here will not reproduce
    d = cpv_periodsearch(
        time, flat_flux, starid, cachedir, t0='binmin', periodogram_method='pdm'
    )

    cd = {
        'method': 'psplinefilt_findpeaks',
        'height': '2_P2P',
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

    outpath = join(cachedir, f'{starid}_phase.png')
    titlestr = starid
    plot_phased_light_curve(
        d['times'], d['fluxs'], d['t0'], d['period'], outpath,
        titlestr=titlestr, binsize_phase=0.005, showtimeticks=None,
        savethefigure=True
    )

if __name__ == "__main__":
    epss = np.arange(0, 0.01 + 1e-3, 1e-3).round(3)
    for eps in epss:
        main(eps)
