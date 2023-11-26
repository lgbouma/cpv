"""
HIRES order 10 showed some weird continuum normalization stuff
"""
# rj order 10...
import pickle
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from numpy import array as nparr
from os.path import join

from scipy.ndimage import gaussian_filter1d
from rudolf.plotting import multiline

def make_plot(
    do_5880 = 1,
    orderstr = 'rj_order10',
    ylim = None
):

    cachedir = '/Users/luke/Dropbox/proj/cpv/results/HIRES_continuum_evoln/'
    # made by plot_lineevolnpanel.py
    if orderstr == 'rj_order10' and do_5880:
        pklpath = join(
            cachedir, 'spec_cache_rj_order10_5880norm.pkl'
        )
    else:
        pklpath = join(
            cachedir, f'spec_cache_{orderstr}.pkl'
        )

    with open(pklpath, 'rb') as f:
        d = pickle.load(f)

    wav = d['wav']
    flxs = d['flxs']
    wavs = [wav for ix in range(len(flxs))]
    mjds = d['mjds']
    flx_median = d['flx_median']

    fn = lambda x: gaussian_filter1d(x, sigma=5)
    flxs = [fn(f) for f in flxs]

    from aesthetic.plot import set_style, savefig
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
    if do_5880:
        ax.set_ylabel('Relative flux (5880A norm)')
    else:
        ax.set_ylabel('Relative flux (order med norm)')

    if isinstance(ylim, list):
        ax.set_ylim(ylim)
    else:
        ax.set_ylim([0,3])

    outdir = '/Users/luke/Dropbox/proj/cpv/results/HIRES_continuum_evoln'
    s = ''
    if do_5880:
        s += '5880norm'
    else:
        s += 'ordermednorm'
    if isinstance(ylim, list):
        s += f'_ylim{ylim[0]}-{ylim[1]}'
    outpath = os.path.join(outdir, f'{orderstr}_continuum_evoln_{s}.png')

    savefig(fig, outpath, dpi=400)

if __name__ == "__main__":
    make_plot(do_5880=1, orderstr='rj_order10')
    make_plot(do_5880=0, orderstr='rj_order10')

    orderstrs = [
        'ij_order00', 'rj_order11', 'bj_order13', 'bj_order07', 'bj_order05',
        'bj_order04'
    ]
    for o in orderstrs:
        make_plot(do_5880=0, orderstr=o)
        make_plot(do_5880=0, orderstr=o, ylim=[0,10])
