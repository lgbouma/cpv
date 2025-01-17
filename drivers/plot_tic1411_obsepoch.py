import pickle
import numpy as np, matplotlib.pyplot as plt, pandas as pd
import os
from os.path import join

from matplotlib import rcParams
from aesthetic.plot import set_style, savefig

# bringing out lots of other projects...
from complexrotators.plotting import (
    prepare_local_lc, plot_phased_light_curve
)
from gyrojo.plotting import add_gradient_patch
from cdips.lcproc import detrend as dtr

outdir = '/Users/luke/Dropbox/proj/cpv/results/tic1411_obsepoch'

ticid = "141146667"

lcpath = '/Users/luke/Dropbox/proj/cpv/data/photometry/tess/mastDownload/TESS/tess2024030031500-s0075-0000000141146667-0270-s/tess2024030031500-s0075-0000000141146667-0270-s_lc.fits'


(time, flux, qual, x_obs, y_obs, y_flat,
 y_trend, x_trend, cadence_sec, sector,
 starid) = prepare_local_lc(ticid, lcpath, None, fluxkey='SAP_FLUX')

sel = ~( (time > 3359.4) & (time < 3360.15) )
time, flux = time[sel], flux[sel]

flux /= np.nanmedian(flux)

# light detrending...
y_flat, y_trend = dtr.detrend_flux(
    time, flux, method='biweight', cval=2, window_length=4.0,
    break_tolerance=0.5
)

flux = 100 * (y_flat - np.nanmedian(y_flat))


#from copy import deepcopy
#time = deepcopy(x_obs)
#flux = deepcopy(y_obs) # lightly detrended


# time and flux are raw PDCSAP here...

# binned to 10 minutes
from astrobase.lcmath import time_bin_magseries
bd = time_bin_magseries(time, flux, binsize=600, minbinelems=1)

# make the plot
plt.close("all")
set_style('science')
rcParams['font.family'] = 'Arial'

fig = plt.figure(figsize=(6,4))
axd = fig.subplot_mosaic(
    """
    AAAA
    AAAA
    BBCC
    BBCC
    BBCC
    """
)
fig = plt.figure(figsize=(6,5))
axd = fig.subplot_mosaic(
    """
    AAAA
    AAAA
    BBCC
    BBCC
    BBCC
    DDDD
    DDDD
    """
)


# the light curve
ax = axd['A']

ax.scatter(time, flux,
           c='lightgray', s=1, zorder=1, rasterized=True)
ax.scatter(bd['binnedtimes'], bd['binnedmags'],
           c='k', s=0.2, zorder=2, rasterized=True)

xmin = 3357.949
xmax = xmin + 5.2/24
ymin, ymax = ax.get_ylim()
add_gradient_patch(ax, xmin, xmax, ymin, ymax, resolution=100,
                   color0='lime', color1='lime')
ax.update({
    'xlim': [3351, 3363.7],
    'ylim': [-10, 6],
    'xlabel': 'Time [BTJD]',
    'ylabel': '$\Delta$ Flux [%]'
})

# before
ax = axd['B']

mask = (time > 3352) & (time < 3356)
period = 0.163762133#*24
t0 = 3339.9326
t0 += 0.15*period
binsize_phase = 0.005
xlim = [-0.1, 1.3]
c0 = 'darkgray'
c1 = 'k'
phasewrap, longwrap = False, True#False

plot_phased_light_curve(
    time[mask], flux[mask]/100, t0, period, None, fig=fig, ax=ax,
    titlestr='BTJD 3352 - 3356', binsize_phase=binsize_phase,
    xlim=xlim, yoffset=0,
    showtext=None,
    savethefigure=False, dy=0, rasterized=True, c0=c0, c1=c1,
    phasewrap=phasewrap, longwrap=longwrap
)
ax.set_ylabel(r"$\Delta$ Flux [%]")
ax.set_xlabel(r"Phase, φ")
ax.set_ylim([-10, 6])
#ax.set_xlim([-0.1, 1.3])

# after
ax = axd['C']

mask = (
    (time > 3358) & (time < 3363)
    & ~( (time > 3359.4) & (time < 3359.6) )
)
plot_phased_light_curve(
    time[mask], flux[mask]/100, t0, period, None, fig=fig, ax=ax,
    titlestr='BTJD 3358 - 3363', binsize_phase=binsize_phase,
    xlim=xlim, yoffset=0, showtext=None,
    savethefigure=False, dy=0, rasterized=True, c0=c0, c1=c1,
    phasewrap=phasewrap, longwrap=longwrap
)
ax.set_ylabel(r"$\Delta$ Flux [%]")
ax.set_xlabel(r"Phase, φ")
ax.set_ylim([-10, 6])
#ax.set_xlim([-0.1, 1.3])

# way after; evolving scatter
ax = axd['D']

ax.scatter(time, flux,
           c='lightgray', s=1, zorder=1, rasterized=True)
ax.scatter(bd['binnedtimes'], bd['binnedmags'],
           c='k', s=0.2, zorder=2, rasterized=True)
ax.update({
    'xlim': [3360.3, 3363.7],
    'ylim': [-10, 6],
    'xlabel': 'Time [BTJD]',
    'ylabel': '$\Delta$ Flux [%]'
})


fig.tight_layout()

outpath = join(outdir, f'tic1411_obsepoch.pdf')
savefig(fig, outpath, dpi=400)

