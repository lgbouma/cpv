from astroARIADNE.star import Star
from astroARIADNE.plotter import SEDPlotter
import os
from os.path import join

starname = 'NGTS-6'
out_folder = f'../results/ariadne_sed_fitting/{starname}'
in_file = os.path.join(out_folder, 'BMA.pkl')
plots_out_folder = join(out_folder, 'plots')
if not os.path.exists(plots_out_folder): os.mkdir(plots_out_folder)

artist = SEDPlotter(in_file, plots_out_folder)
artist.plot_SED_no_model()
artist.plot_SED()
artist.plot_bma_hist()
artist.plot_bma_HR(10)
artist.plot_corner()


