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
from matplotlib.transforms import blended_transform_factory

from astropy.io import fits

from scipy.ndimage import gaussian_filter1d

from rudolf.plotting import multiline
from cdips_followup.spectools import read_hires, LINE_D
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
        try:
            flx_2d, wav_2d = read_hires(
                specpath, is_registered=0, return_err=0, start=10, end=-10
            )
        except IndexError:
            # odd edge case for j537.174 blue chip...
            print(f'caught index error and pulling just flux for {specpath}...')
            hdul = fits.open(specpath)
            flx_2d = hdul[0].data
            flx_2d = flx_2d[:, 10:-10]
            hdul.close()

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

    # TIC 4029 (j533 also exists, and is decent)
    datadir = "/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC402980664_RDX"
    outdir = '/Users/luke/Dropbox/proj/cpv/results/HIRES_results'
    datestr = 'j547'

    #  # TIC 1411
    #  datadir = "/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/TIC141146667_RDX"
    #  outdir = '/Users/luke/Dropbox/proj/cpv/results/HIRES_results'
    #  datestr = 'j537'

    chips = 'b,r,i'.split(",")
    norders = [23,16,10]

    for chip,norder in zip(chips,norders):

        spectrum_paths = np.sort(glob(join(datadir, f'{chip}{datestr}*fits')))

        if 'TIC1411' in datadir:
            # trim few junky outliers
            spectrum_paths = spectrum_paths[1:-2]
        if 'TIC4029' in datadir and datestr == 'j533':
            spectrum_paths = spectrum_paths[:-1]
        if 'TIC4029' in datadir and datestr == 'j547':
            spectrum_paths = spectrum_paths[1:]

        for order in range(norder):

            orderstr = f"{chip}j_order{str(order).zfill(2)}"

            mjds, mednorm_flxs, wavs = given_specpaths_get_timeseries_data(
                spectrum_paths, order
            )

            fn = lambda x: gaussian_filter1d(x, sigma=5)
            flxs = [fn(f) for f in mednorm_flxs]

            # make plot
            plt.close("all")
            set_style("science")

            fig, ax = plt.subplots(figsize=(10,3))
            lc = multiline(
                wavs, flxs, 24*(nparr(mjds)-np.min(mjds)),
                cmap='viridis',
                ax=ax, lw=0.5
            )

            loc = 'lower left' if not isinstance(ylim, (list, tuple)) else 'upper right'
            axins1 = inset_axes(ax, width="20%", height="5%", loc=loc, borderpad=1.5)
            cb = fig.colorbar(lc, cax=axins1, orientation="horizontal")
            cb.ax.tick_params(labelsize='xx-small')
            cb.ax.set_title('Time [hours]', fontsize='xx-small')
            cb.ax.tick_params(size=0, which='both') # remove the ticks
            axins1.xaxis.set_ticks_position("bottom")

            ax.set_xlabel('λ [$\AA$]')
            ax.set_ylabel('Relative flux (order med norm)')

            if isinstance(ylim, (list,tuple)):
                ax.set_ylim(ylim)
            else:
                ax.set_ylim([0,3])

            # plot line names
            xmin = min(wavs[0])
            xmax = max(wavs[0])
            this_d = []
            for k, v in LINE_D:
                if v > xmin and v<xmax:
                    this_d.append([k, v])
            if len(this_d) > 0:
                for k, v in this_d:
                    ylim = ax.get_ylim()
                    delta_y = 0.9*(max(ylim) - min(ylim))
                    ax.vlines(v, min(ylim)+delta_y, max(ylim), zorder=-3,
                              linestyles=':', color='k', lw=0.3)
                    ax.set_ylim(ylim)

                    tform = blended_transform_factory(ax.transData, ax.transAxes)
                    ax.text(v, 0.95, k, ha='center', va='top', transform=tform,
                            fontsize=4)


            outdir = f'/Users/luke/Dropbox/proj/cpv/results/HIRES_continuum_evoln/{datestr}'
            if not os.path.exists(outdir): os.mkdir(outdir)
            s = ''
            s += 'ordermednorm'
            if isinstance(ylim, (list,tuple)):
                s += f'_ylim{ylim[0]}-{ylim[1]}'
            outpath = os.path.join(outdir, f'{orderstr}_continuum_evoln_{s}.png')

            #plt.show(block=True)

            savefig(fig, outpath, dpi=400, writepdf=0)

if __name__ == "__main__":

    plot_continuum_timeseries(ylim=list([0,10]))
    plot_continuum_timeseries(ylim=list([0,3]))
