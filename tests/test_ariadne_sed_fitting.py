"""
Demo for astroARIADNE
https://github.com/jvines/astroARIADNE
Vines & Jenkins 2022

environment: (py38_ariadne)
"""

##################
# query the star #
##################
from astroARIADNE.star import Star
import os
from os.path import join

ra = 75.795
dec = -30.399
starname = 'NGTS-6'
gaia_id = 4875693023844840448

s = Star(starname, ra, dec, g_id=gaia_id)

######################
# fit the photometry #
######################
from astroARIADNE.fitter import Fitter

out_folder = f'../results/ariadne_sed_fitting/{starname}'

engine = 'dynesty'
nlive = 500
dlogz = 0.5
bound = 'multi'
sample = 'rwalk'
threads = 4
dynamic = False

setup = [engine, nlive, dlogz, bound, sample, threads, dynamic]

# Feel free to uncomment any unneeded/unwanted models
models = [
	#'phoenix',
	#'btsettl',
	#'btnextgen',
	#'btcond',
	'kurucz',
	#'ck04'
]

f = Fitter()
f.star = s
f.setup = setup
f.av_law = 'fitzpatrick'
f.out_folder = out_folder
f.bma = True
f.models = models
f.n_samples = 100000

# The default prior for Teff is an empirical prior drawn from the RAVE survey
# temperatures distribution, the distance prior is drawn from the Bailer-Jones
# distance estimate from Gaia EDR3, and the radius has a flat prior ranging from
# 0.5 to 20 R$_\odot$. The default prior for the metallicity z and log g are also
# their respective distributions from the RAVE survey, the default prior for Av
# is a flat prior that ranges from 0 to the maximum of line-of-sight as per the
# SFD map, finally the excess noise parameters all have gaussian priors centered
# around their respective uncertainties.
f.prior_setup = {
        'teff': ('default'),
        'logg': ('default'),
        'z': ('default'),
        'dist': ('default'),
        'rad': ('default'),
        'Av': ('default')
}

# So if you knew (from a spectroscopic analysis, for example) that the effective
# temperature is 5600 +/- 100 and the metallicity is [Fe/H] = 0.09 +/- 0.05 and
# you wanted to use them as priors, and the star is nearby (< 70 pc), so you
# wanted to fix Av to 0, your prior dictionary should look like this:
f.prior_setup = {
        'teff': ('normal', 4800, 300),
        'logg': ('default'),
        'z': ('normal', 0.1, 0.1),
        'dist': ('default'),
        'rad': ('default'),
        'Av': ('fixed', 0)
}

f.initialize()
# this takes like 10 minutes
f.fit_bma()

##############
# make plots #
##############
starname = 'NGTS-6'
out_folder = f'../results/ariadne_sed_fitting/{starname}'
in_file = os.path.join(out_folder, 'BMA.pkl')
plots_out_folder = join(out_folder, 'plots')
if not os.path.exists(plots_out_folder): os.mkdir(plots_out_folder)

from astroARIADNE.plotter import SEDPlotter
artist = SEDPlotter(in_file, plots_out_folder)
artist.plot_SED_no_model()
artist.plot_SED()
artist.plot_bma_hist()
artist.plot_bma_HR(10)
artist.plot_corner()

import IPython; IPython.embed()
