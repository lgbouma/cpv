from glob import glob
from os.path import join
import os
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from cdips_followup.spectools import viz_1d_spectrum
from cdips_followup.spectools import read_hires

datadir = "/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC141146667_RDX"
specpaths = glob(join(datadir, "ij*.fits"))
idstring = 'TIC141146667'
outdir = '/Users/luke/Dropbox/proj/cpv/results/TIC141146667_HIRES_raw'
if not os.path.exists(outdir): os.mkdir(outdir)
datestr = 'j537'

for specpath in np.sort(specpaths):

    flx_2d, wav_2d = read_hires(specpath, is_registered=0, return_err=0)

    start = 10
    end = -10
    order = 0 # halpha
    flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]

    outname = f'{idstring}_{os.path.basename(specpath).rstrip(".fits")}_order{str(order).zfill(2)}.png'
    outpath = join(outdir, outname)

    viz_1d_spectrum(flx, wav, outpath, xlim=[6542,6582], vlines=None, names=None,
                    ylim=[200,1300], norm_median=False, ylabel=None, xlabel=None,
                    fig=None, ax=None, axtitle=None)

    wav0 = 6562.8
    outname = f'{idstring}_{os.path.basename(specpath).rstrip(".fits")}_order{str(order).zfill(2)}_showvel.png'
    outpath = join(outdir, outname)
    viz_1d_spectrum(flx, wav, outpath, xlim=[6542,6582], vlines=None, names=None,
                    ylim=[200,1300], norm_median=False, ylabel=None, xlabel=None,
                    fig=None, ax=None, axtitle=None, show_vel=1, wav0=wav0)
