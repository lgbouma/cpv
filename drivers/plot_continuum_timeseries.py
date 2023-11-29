"""
Given a set of time-series spectra, plot all the orders, over time.
"""
import pickle
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from numpy import array as nparr
from os.path import join
from glob import glob

from astropy.io import fits

from scipy.ndimage import gaussian_filter1d

from rudolf.plotting import multiline
from cdips_followup.spectools import read_hires
from aesthetic.plot import set_style, savefig

def given_specpaths_get_timeseries_data(spectrum_paths, order):

    mjds, mednorm_flxs, wavs = [],[],[]
    for ix, specpath in enumerate(spectrum_paths):

        # header time info
        hdul = fits.open(specpath)
        dateobs = hdul[0].header['DATE-OBS']
        exptime = hdul[0].header['EXPTIME']
        mjd = hdul[0].header['MJD']
        utc = hdul[0].header['UTC'][:5]
        timestr = f"{dateobs} {utc}UT"
        hdul.close()

        # flux info
        flx_2d, wav_2d = read_hires(
            specpath, is_registered=0, return_err=0, start=10, end=-10
        )

        flx, wav = flx_2d[order, :], wav_2d[order, :]

        # median normalize over order
        flx /= np.nanmedian(flx)

        mjds.append(mjd)
        mednorm_flxs.append(flx)
        wavs.append(wav)

    mjds = np.array(mjds).astype(float)
    mednorm_flxs = np.array(mednorm_flxs)
    wavs = np.array(wavs)

    # does the wavelength solution vary over time?  it should not.  this
    # assertion statement verifies that.
    assert np.diff(wavs, axis=0).sum() == 0

    return mjds, mednorm_flxs, wavs



def plot_continuum_timeseries(ylim=None):

    datadir = "/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC402980664_RDX"
    outdir = '/Users/luke/Dropbox/proj/cpv/results/HIRES_results'

    chips = 'b,r,i'.split(",")
    norders = [23,16,10]

    for chip,norder in zip(chips,norders):

        spectrum_paths = np.sort(glob(join(datadir, f'{chip}j*fits')))

        for order in range(norder):

            orderstr = f"{chip}j_order{str(order).zfill(2)}"

            mjds, mednorm_flxs, wavs = given_specpaths_get_timeseries_data(
                spectrum_paths, order
            )

            fn = lambda x: gaussian_filter1d(x, sigma=5)
            flxs = [fn(f) for f in mednorm_flxs]

            # make plot
            plt.close("all")
            set_style("clean")

            fig, ax = plt.subplots(figsize=(10,3))
            lc = multiline(
                wavs, flxs, 24*(nparr(mjds)-np.min(mjds)),
                cmap='viridis',
                ax=ax, lw=0.5
            )

            loc = 'lower left' if not isinstance(ylim, list) else 'upper left'
            axins1 = inset_axes(ax, width="20%", height="5%", loc=loc, borderpad=1.5)
            cb = fig.colorbar(lc, cax=axins1, orientation="horizontal")
            cb.ax.tick_params(labelsize='xx-small')
            cb.ax.set_title('Time [hours]', fontsize='xx-small')
            cb.ax.tick_params(size=0, which='both') # remove the ticks
            axins1.xaxis.set_ticks_position("bottom")

            ax.set_xlabel('λ [$\AA$]')
            ax.set_ylabel('Relative flux (order med norm)')

            if isinstance(ylim, list):
                ax.set_ylim(ylim)
            else:
                ax.set_ylim([0,3])

            outdir = '/Users/luke/Dropbox/proj/cpv/results/HIRES_continuum_evoln'
            s = ''
            s += 'ordermednorm'
            if isinstance(ylim, list):
                s += f'_ylim{ylim[0]}-{ylim[1]}'
            outpath = os.path.join(outdir, f'{orderstr}_continuum_evoln_{s}.png')

            savefig(fig, outpath, dpi=400)

if __name__ == "__main__":

    plot_continuum_timeseries()
