"""
Given the TIC1411 HIRES time-series data, fit a model for the emission
components.

Given the complexity of the data, it is challenging to know how complex of a
model to specify, or what free parameters it should be allowed to have.

The emission at v/v_eq > looks quite sinusoidal.  However, the amplitude of the
emission is time-variable.  It also not clear whether to separate the two
components, or not.

So, perhaps consider models where the emission in flux is a sum of gaussians (or a
lorentzians, or voigts), with means

    μ_0 ~ K_0(t) sin(ω0t + φ0)
    μ_1 ~ K_1(t) sin(ω1t + φ1)

and time-varying amplitudes A_n(t).
"""

import os, pickle
from os.path import join
from glob import glob
from datetime import datetime
from copy import deepcopy
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from numpy import array as nparr

from complexrotators.paths import RESULTSDIR
plotdir = join(RESULTSDIR, "emission_modelfitter")
if not os.path.exists(plotdir): os.mkdir(plotdir)

from complexrotators.modelfitter import ModelFitter

# Pre-requisite: running plot_movie_specriver.py
# NOTE: if you change v_eq, this requires a rerun.
pklpath = (
    '/Users/luke/Dropbox/proj/cpv/results/movie_specriver/'
    'TIC_141146667/141146667_specriver_cache.pkl'
)
with open(pklpath, 'rb') as f:
    d = pickle.load(f)

yval = 1.*d['xvals'][0]
xval = (d['spectimes']-d['t0'])/d['period']
xval -= np.min(np.ceil(xval))
flux_arr = d['spec_flux_arr']

mask = np.abs(yval) < 1.2

flux_arr_masked = deepcopy(flux_arr) - 1
flux_arr_maskedinv = deepcopy(flux_arr) -1

flux_arr_masked[mask, :] = np.nan  # Set masked values to NaN
flux_arr_maskedinv[~mask, :] = np.nan  # Set non-masked values to NaN

# make err array based on point-to-point variation outside line core
from cdips.utils.lcutils import p2p_rms
err_flux_arr = np.zeros_like(flux_arr)
for ix in range(flux_arr.shape[1]):
    err_flux_arr[:,ix] = p2p_rms(flux_arr_masked[:,ix])


# viz
plt.close("all")
fig, ax = plt.subplots()
cmap = plt.get_cmap('Greys')
cmap.set_bad(color='white')  # Set color for NaN values
c = ax.pcolormesh(xval, yval, flux_arr_masked, cmap=cmap, shading='auto', rasterized=True)
cb = fig.colorbar(c, ax=ax)
cb.set_label('f (med norm. continuum)')
ax.update({'xlabel':'time (P)', 'ylabel': 'Δv/veq'})
fig.savefig(
    join(plotdir, "basic_med_norm_masked_tic1411_j537.png"), dpi=300,
    bbox_inches='tight'
)

# viz
plt.close("all")
fig, ax = plt.subplots()
cmap = plt.get_cmap('Greys')
cmap.set_bad(color='white')  # Set color for NaN values
c = ax.pcolormesh(xval, yval, flux_arr_maskedinv, cmap=cmap, shading='auto', rasterized=True)
cb = fig.colorbar(c, ax=ax)
cb.set_label('f (med norm. continuum)')
ax.update({'xlabel':'time (P)', 'ylabel': 'Δv/veq'})
fig.savefig(
    join(plotdir, "basic_med_norm_maskedinv_tic1411_j537.png"), dpi=300,
    bbox_inches='tight'
)

# define guessdict
N_gaussians = 3
for N_gaussians in [2,3]:

    if N_gaussians == 2:
        K_n_guess = nparr([2.3, 4])
        phi_n_guess = nparr([0.1, 3.1])
        sigma_n_guess = nparr([0.1, 0.1])
        A_n_guess = nparr([0.4, 0.4])
        period_guess = nparr([1.0]*2)
    elif N_gaussians == 3:
        K_n_guess = nparr([2.1, 2.6, 4])
        phi_n_guess = nparr([0.1, 0.1, 3.1])
        sigma_n_guess = nparr([0.1, 0.1, 0.1])
        A_n_guess = nparr([0.2, 0.2, 0.4])
        period_guess = nparr([1.0]*3)

    guessdict = {
        'K': K_n_guess,
        'phi': phi_n_guess,
        'sigma': sigma_n_guess,
        'A': A_n_guess,
        'period': period_guess,
    }

    # fit data
    m = ModelFitter(
        xval, yval, flux_arr_masked.T, err_flux_arr.T,
        modelid=f'{N_gaussians}_gaussians', N_samples=1000, N_cores=2, N_chains=2,
        plotdir=plotdir, overwrite=False, guessdict=guessdict,
        map_guess_method='handtuned', flux_arr_nomask=flux_arr.T
    )
