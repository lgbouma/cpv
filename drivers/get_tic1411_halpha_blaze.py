"""
The TIC1411 sequence default reduction has an evolving blaze function.  Turn
this off.  Adopt the "median"est blaze function as your actual one for Hα
analyses.
"""
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from astropy.io import fits
from glob import glob
import os
from os.path import join
from cdips_followup.spectools import read_hires

def get_ylimguess(y):
    ylow = np.nanpercentile(y, 2.5)
    yhigh = np.nanpercentile(y, 99.5)
    ydiff = (yhigh-ylow)
    ymin = ylow - 0.15*ydiff
    ymax = yhigh + 0.35*ydiff
    return [ymin,ymax]

deblzdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC141146667_DEBLAZED/data/mir3/iodfitsdb'
reducdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC141146667_REDUCED/data/mir3/iodfits'

deblzpaths = np.sort(glob(join(deblzdir, 'ij*fits')))
reducpaths = np.sort(glob(join(reducdir, 'ij*fits')))

specpath = '/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC141146667_DEBLAZED/data/mir3/iodfitsdb/ij537.172.fits'
_, wav_2d = read_hires(specpath, is_registered=0, return_err=0)

blz_flxs = []
idxs = []

for deblzpath, reducpath in zip(deblzpaths, reducpaths):

    idx = os.path.basename(deblzpath)

    try:
        deblz_flx_2d, wav_2d = read_hires(deblzpath, is_registered=0, return_err=0)
    except IndexError:
        # odd edge case for j537.174 blue chip...
        print('caught index error and pulling just flux...')
        hdul = fits.open(deblzpath)
        deblz_flx_2d = hdul[0].data
        hdul.close()

    hdul = fits.open(reducpath)
    reduc_flx_2d = hdul[0].data
    hdul.close()

    blz_flx_2d = reduc_flx_2d / deblz_flx_2d

    start = 10
    end = -10
    order = 0
    blz_flx, wav = blz_flx_2d[order, start:end], wav_2d[order, start:end]
    reduc_flx = reduc_flx_2d[order, start:end]
    deblz_flx = deblz_flx_2d[order, start:end]

    blz_flxs.append(blz_flx)
    idxs.append(idx)

    outdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC141146667_BLAZE'
    savpath = join(
        outdir, os.path.basename(deblzpath).replace(".fits","_viz.png")
    )
    if not os.path.exists(savpath):
        plt.close('all')
        fig, axs = plt.subplots(nrows=3, figsize=(12,6), sharex=True)
        axs[0].plot(wav, deblz_flx, c='k', lw=0.5)
        axs[0].set_ylabel('deblz_flx')
        axs[0].set_ylim(get_ylimguess(deblz_flx))
        axs[1].plot(wav, reduc_flx, c='k', lw=0.5)
        axs[1].set_ylabel('reduc_flx')
        axs[1].set_ylim(get_ylimguess(reduc_flx))
        axs[2].plot(wav, blz_flx, c='k', lw=0.5)
        axs[2].set_ylabel('blz')
        axs[2].set_ylim((0.5, 1.4))
        fig.tight_layout()
        fig.savefig(savpath, bbox_inches='tight', dpi=300)
        print(f'wrote {savpath}')

plt.close("all")

fig, ax = plt.subplots(figsize=(12,2.5))

blz_flxs = np.array(blz_flxs)
med_blz_flx = np.nanmedian(blz_flxs, axis=0)
differences = np.abs(blz_flxs - med_blz_flx)
total_differences = np.nansum(differences, axis=1)
most_median_like_sample_index = np.nanargmin(total_differences)

for idx, blz_flx in zip(idxs, blz_flxs):
    if idx not in [ 'ij537.169.fits', 'ij537.165.fits' ] and '537' in idx:
        ax.plot(wav, blz_flx, c='k', alpha=0.2, lw=0.5)

ax.plot(wav, med_blz_flx, c='C0', lw=0.5)
ax.plot(wav, blz_flxs[most_median_like_sample_index,:], c='C1', lw=1)

outdf = pd.DataFrame({
    'wav': wav,
    'blz_flx': blz_flxs[most_median_like_sample_index,:]
})
outdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC141146667_BLAZE'

csvpath = join(
    outdir, "j537_ichip_order00_preferred_median_blaze_flux.csv"
)
outdf.to_csv(csvpath, index=False)

savpath = join(
    outdir, "j537_ichip_order00_median_blaze_flux.png"
)
fig.savefig(savpath, bbox_inches='tight', dpi=300)

