import numpy as np, matplotlib.pyplot as plt
from os.path import join
import batman

from complexrotators.lcprocessing import (
    cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
)
from complexrotators.plotting import (
    plot_phased_light_curve, plot_dipcountercheck, plot_cpvvetter
)


P0 = 0.321

def get_ytra_ysec(t):

    params = batman.TransitParams()
    params.t0 = 0.                       #time of inferior conjunction
    params.per = P0                      #orbital period
    params.rp = 0.2                      #planet radius (in units of stellar radii)
    params.a = 2                       #semi-major axis (in units of stellar
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
    A0 = 0.05
    phi1 = 1.12
    A1 = 0.05

    y0 = A0 * np.sin(2*np.pi*t/P0 + phi0)
    y1 = A1 * np.sin(2*np.pi*t/P1 + phi1)

    return y0, y1


def main():

    #
    # generate fake adversarial example
    #
    t = np.linspace(0,10, 5000) # underlying model
    tobs = np.arange(0,10, 2/(60*24)) # tess 120-second cadence

    ytra, ysec = get_ytra_ysec(t)
    y0, y1 = get_y0_y1(t)
    ytot = 1 + y0 + ytra + ysec + y1

    ytra_obs, ysec_obs = get_ytra_ysec(tobs)
    y0_obs, y1_obs = get_y0_y1(tobs)

    np.random.seed(42)
    yerr = np.random.normal(loc=0, scale=0.02, size=len(tobs))

    yobs = 1 + y0_obs + ytra_obs + ysec_obs + y1_obs + yerr

    #
    # plot the truth, and the observed fake data
    #
    fig, axs = plt.subplots(nrows=3, figsize=(10,6))

    axs[0].plot(t, ytot)
    axs[1].scatter(tobs, yobs, c='k', s=1, linewidths=0)
    ys = [y0+ytra+ysec, y1, ytra+ysec]
    labels = ['P0: sinusoid plus eclipses', 'P1: 3/2*P0', 'just eclipses at P0']
    for ix, (y, l) in enumerate(zip(ys, labels)):
        axs[2].scatter(t, y, c=f'C{ix}', s=1, linewidths=0, label=l)
    axs[2].legend(loc='best', fontsize='small')
    axs[2].set_xlabel('time [days]')

    fig.savefig('testadversarial321.png', dpi=400, bbox_inches='tight')
    plt.close("all")

    #
    # process fake observed data to see if we would be fooled
    #
    starid = "testadversarial321"
    cachedir = "./"

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
    main()
