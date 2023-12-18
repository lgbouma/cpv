import os
import pandas as pd, matplotlib.pyplot as plt, numpy as np
from glob import glob
from os.path import join
from complexrotators.paths import RESULTSDIR

from aesthetic.plot import set_style, savefig

def t_to_phase(t, dofloor=0):
    t0 = 2450000 + 1791.12 + 0.5*18.5611/24
    period = 18.5611/24
    φ = (t-t0)/period
    if dofloor:
        φ = (t-t0)/period- np.floor( (t-t0)/period )
    return φ

def get_groundphot(
    inst='KeplerCam',
    bandpass='g',
    utcdatestr='20231208'
):

    assert utcdatestr in ['20231208','20231215','20231216']
    assert inst in ['KeplerCam','Tierras']
    assert bandpass in [None, 'g', 'B', 'Tierras']

    TIERRASDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/TIERRAS'
    KEPLERCAMDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/KeplerCam'

    if inst == 'Tierras':
        csvpath = glob(join(
            TIERRASDIR,
            f"{utcdatestr}_TIC402980664_circular_fixed_ap_phot_*.csv")
        )[0]
        df = pd.read_csv(csvpath)
        df['Target Relative Flux'] *= 10

    elif inst == 'KeplerCam':

        df = pd.read_csv(
            join(KEPLERCAMDIR, f"LP12-502_{utcdatestr}_KeplerCam_{bandpass}.dat"),
            delim_whitespace=True
        )
        df['BJD TDB'] = df['BJD_TDB_B']
        df['Target Relative Flux'] = df['rel_flux_T1_n']
        df['Target Relative Flux Error'] = df['rel_flux_err_T1_n']

    else:
        raise NotImplementedError

    t = np.array(df['BJD TDB'])
    phase = t_to_phase(t, dofloor=1)
    flux = np.array(df['Target Relative Flux'])
    flux_err = np.array(df['Target Relative Flux Error'])

    return t, phase, flux, flux_err

def main(insts, bandpasses, utcdatestrs):

    bandcolors = {
        'B': 'C3',
        'g': 'C2',
        'Tierras': 'C1',
    }
    datemarkers = {
        '20231208': 'o',
        '20231215': 'X',
        '20231216': 's'
    }

    set_style('clean')
    fig, ax = plt.subplots(figsize=(4,3))

    for i,b,u in zip(insts, bandpasses, utcdatestrs):

        t, phase, flux, flux_err = get_groundphot(
            inst=i, bandpass=b, utcdatestr=u
        )

        yval = flux - np.nanmedian(flux)
        yerr = flux_err

        l = f"{u}: {i} {b}"

        ax.scatter(phase, 1e2*yval, c=bandcolors[b], alpha=0.9,
                   linewidths=0, zorder=2, s=2, marker=datemarkers[u], label=l)
        ax.errorbar(phase, 1e2*yval, yerr=1e2*yerr,
                    c=bandcolors[b], ecolor=bandcolors[b], elinewidth=0.5,
                    lw=0, alpha=0.1, zorder=1)

        ax.set_ylabel('$f$ [%]')
        ax.set_ylim([-1e2*0.07, 1e2*0.04])

    ax.legend(loc='lower right', fontsize='xx-small')
    ax.set_xlabel("Phase, φ")

    outdir = join(RESULTSDIR, 'chromaticity_keplercam_tierras')
    if not os.path.exists(outdir): os.mkdir(outdir)

    uinsts = "_".join(np.unique(insts))
    udates = "_".join(np.unique(utcdatestrs))
    outpath = join(outdir, f'flux_vs_phase_{uinsts}_{udates}.png')
    savefig(fig, outpath)

if __name__ == "__main__":

    insts = ['KeplerCam', 'KeplerCam', 'KeplerCam', 'KeplerCam', 'Tierras',
             'KeplerCam', 'KeplerCam']
    bandpasses = ['B', 'g', 'B', 'g', 'Tierras', 'B', 'g']
    utcdatestrs = ['20231208', '20231208', '20231215', '20231215', '20231215',
                   '20231216', '20231216']
    main(insts, bandpasses, utcdatestrs)

    insts = ['KeplerCam', 'KeplerCam']
    bandpasses = ['B', 'g']
    utcdatestrs = ['20231208', '20231208']
    main(insts, bandpasses, utcdatestrs)

    insts = ['KeplerCam', 'KeplerCam', 'Tierras']
    bandpasses = ['B', 'g', 'Tierras']
    utcdatestrs = ['20231215', '20231215', '20231215']
    main(insts, bandpasses, utcdatestrs)

    insts = ['KeplerCam', 'KeplerCam']
    bandpasses = ['B', 'g']
    utcdatestrs = ['20231216', '20231216']
    main(insts, bandpasses, utcdatestrs)

    insts = ['KeplerCam', 'KeplerCam', 'KeplerCam', 'KeplerCam', 'Tierras']
    bandpasses = ['B', 'g', 'B', 'g', 'Tierras']
    utcdatestrs = ['20231208', '20231208', '20231215', '20231215', '20231215']
    main(insts, bandpasses, utcdatestrs)

