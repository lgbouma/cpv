import os, pickle
from os.path import join
from glob import glob
from datetime import datetime
from copy import deepcopy
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from numpy import array as nparr

from complexrotators.paths import RESULTSDIR
plotdir = join(RESULTSDIR, "test_emission_modelfitter")
if not os.path.exists(plotdir): os.mkdir(plotdir)

from complexrotators.modelfitter import (
    ModelFitter, generate_synthetic_data_specific_params
)

# Define specific true parameters
num_time_points = 50
K_n_true = [2.1, 2.6, 4]
phi_n_true = [0.1, 0.1, 3.1]
sigma_n_true = [0.2, 0.2, 0.2]
A_n_true = [0.2, 0.2, 0.1]  # Constant amplitude
period_true = 1.0
period_true_arr = [1.0]*3

# Generate synthetic data
xval, yval, flux_arr, err_flux_arr = generate_synthetic_data_specific_params(
    K_n_true, phi_n_true, sigma_n_true, A_n_true,
    num_time_points,
    period_true=period_true,
    num_velocity_points=100,
    noise_level=0.1,
    random_seed=42
)

mask = np.abs(yval) < 1.
flux_arr_masked = deepcopy(flux_arr)
flux_arr_masked[:,mask] = np.nan  # Set masked values to NaN

# fit dataaaa

np.random.seed(42)
eps = 0.01
K_n_guess = nparr(K_n_true) + np.random.normal(0, eps)
phi_n_guess = nparr(phi_n_true) + np.random.normal(0, eps)
sigma_n_guess = nparr(sigma_n_true) + np.random.normal(0, eps)
A_n_guess = nparr(A_n_true) + np.random.normal(0, eps)  # Constant amplitude
period_guess = nparr(period_true_arr) + np.random.normal(0, eps)
guessdict = {
    'K': K_n_guess,
    'phi': phi_n_guess,
    'sigma': sigma_n_guess,
    'A': A_n_guess,
    'period': period_guess,
}

m = ModelFitter(
    xval, yval, flux_arr_masked, err_flux_arr,
    modelid='3_gaussians', N_samples=1000, N_cores=2, N_chains=2,
    plotdir=plotdir, overwrite=False, guessdict=guessdict, map_guess_method='handtuned'
)

