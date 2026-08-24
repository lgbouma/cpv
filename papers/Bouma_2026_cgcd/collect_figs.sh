r=$HOME"/Dropbox/proj/cpv/results"
paperdir=`pwd -P`

###############################################

# f1: the CPVs
cp $r/lc_mosaic/lc_mosaic_mosaic4.pdf lc_mosaic_mosaic4.pdf

# fN: MagE mosaic
cp $r/mage_mosaic_20260820/mage_mosaic_simult_tworow_fournightavg.pdf mage_mosaic_simult_tworow_fournightavg.pdf

# fN: MIKE
cp $r/mike_mosaic_20260820/mike_mosaic.pdf mike_mosaic.pdf 

# fN: DBSP mosaic
cp $r/EW_results/DBSP_HAlpha_mosaic_20260818.pdf DBSP_HAlpha_mosaic.pdf

# fN: LP 12-502 DBSP EW profile variability
cp $r/movie_EW_results_showvel_20260818/20231111_DBSP_0016_4.72407671_0.432150_stackew_vs_phase_DBSP_20231111_20231112_showvel.pdf dbsp_20231111_20231112_frame.pdf
blogdir=$HOME"/Dropbox/Documents/websites/portfolio/lgbouma_professional/static/movies"
cp $r/movie_EW_results_showvel_20260818/movie_EW_showvel_20231111_20231112.mov $blogdir/.
cd $blogdir
./convert_mov2mp4.sh
rm movie_EW_showvel_20231111_20231112.mov
mv $blogdir/movie_EW_showvel_20231111_20231112.mp4 $blogdir/LP_12-502_dbsp_Nov2023.mp4
cd $paperdir
cp $blogdir/LP_12-502_dbsp_Nov2023.mp4 figure_additions/.

# fN: FourStar vs TESS
cp $r/fourstar_vs_tess/fourstar_vs_tess_phase.pdf fourstar_vs_tess_phase.pdf

# fN: MuSCAT vs TESS
cp $r/multicolor_phase/TIC262400835_multicolor_phase_stacked_reextracted.pdf TIC262400835_multicolor_phase_stacked.pdf

# fN: KeplerCam & Tierras vs TESS
cp $r/multicolor_phase/LP12502_multicolor_phase_stacked_Dec08Dec15.pdf LP12502_multicolor_phase_stacked_Dec08Dec15.pdf

# fN: empirical power law fits to dip depth data
cp $r/beta_fit/combined/with_sigma_int/money_plot_beta_with_sigma_int.pdf beta_fits.pdf

# fN: cartoon origin scenarios
cp $r/diluted_transmission_spectrum/diluted_transmission_spectrum.pdf diluted_transmission_spectrum.pdf

# fN: cartoon origin scenarios
cp $r/cgcd_cartoons/cgcd_cartoons_b_mixed_geom.pdf cgcd_cartoons_b_mixed_geom.pdf 

###############################################

# supp

# DBSP no dip nights
cp $r/EW_results/DBSP_HAlpha_mosaic_supp.pdf .
cp $r/EW_results/LP12-502_Nov2024_dbsp_vs_tess.pdf .

# HIRES no dip nights
cp $r/hires_lp12502_mosaic/hires_lp12502_mosaic.pdf .

# Dip fitting mosaics
cp $r/dip_mosaic/dip_fit_mosaic_*pdf .
