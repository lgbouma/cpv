from glob import glob
from astropy.io import fits
from os.path import join
import os
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from cdips_followup.spectools import viz_1d_spectrum
from cdips_followup.spectools import read_hires

datadir = "/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC141146667_RDX"
idstring = 'TIC141146667'
outdir = '/Users/luke/Dropbox/proj/cpv/results/TIC141146667_HIRES_raw'
if not os.path.exists(outdir): os.mkdir(outdir)

linestrs, chips, orders, xlims, ylims, wav0s = (
    ['Hα', 'Hγ', 'CaK'],
    ['i', 'b', 'b'],
    [0, 15, 6],
    [[6542, 6582], [4340.5-20, 4340.5+20], [3933.66-20,3933.66+20]],
    [[200, 1300], [-10, 100], [-10, 90]],
    [6562.8, 4340.47, 3933.66]
)

for line, chip, order, xlim, ylim, wav0 in zip(
    linestrs, chips, orders, xlims, ylims, wav0s
):

    specpaths = glob(join(datadir, f"{chip}j*.fits"))
    datestr = 'j537'

    for specpath in np.sort(specpaths):

        try:
            flx_2d, wav_2d = read_hires(specpath, is_registered=0, return_err=0)
        except IndexError:
            # odd edge case for j537.174 blue chip...
            print('caught index error and pulling just flux...')
            hdul = fits.open(specpath)
            flx_2d = hdul[0].data
            hdul.close()

        start = 10
        end = -10

        flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]

        SHOW_VS_WAV = 0
        if SHOW_VS_WAV:
            outname = (f'{idstring}_{line}_'
                       f'{os.path.basename(specpath).rstrip(".fits")}'
                       f'_order{str(order).zfill(2)}.png')
            outpath = join(outdir, outname)

            viz_1d_spectrum(flx, wav, outpath, xlim=xlim, vlines=None, names=None,
                            ylim=ylim, norm_median=False, ylabel=None, xlabel=None,
                            fig=None, ax=None, axtitle=None)

        outname = (f'{idstring}_{line}_'
                   f'{os.path.basename(specpath).rstrip(".fits")}_'
                   f'order{str(order).zfill(2)}_showvel.png')
        outpath = join(outdir, outname)
        if os.path.exists(outpath):
            continue
        viz_1d_spectrum(flx, wav, outpath, xlim=xlim, vlines=None, names=None,
                        ylim=ylim, norm_median=False, ylabel=None, xlabel=None,
                        fig=None, ax=None, axtitle=None, show_vel=1, wav0=wav0)
