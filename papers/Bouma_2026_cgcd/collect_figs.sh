r=$HOME"/Dropbox/proj/cpv/results"
paperdir=`pwd -P`

###############################################

# f1: the CPVs
cp $r/lc_mosaic/lc_mosaic_mosaic4.pdf f1.pdf

# f2: MagE mosaic
cp $r/mage_mosaic_20260820/mage_mosaic_simult_tworow_fournightavg.pdf f2.pdf

# f3: MIKE
cp $r/mike_mosaic_20260820/mike_mosaic.pdf f3.pdf 

# f4: DBSP mosaic
cp $r/EW_results/DBSP_HAlpha_mosaic_20260818.pdf f4.pdf

# f5: LP 12-502 DBSP EW profile variability
cp $r/movie_EW_results_showvel_20260818/20231111_DBSP_0016_4.72407671_0.432150_stackew_vs_phase_DBSP_20231111_20231112_showvel.pdf f5.pdf
blogdir=$HOME"/Dropbox/Documents/websites/portfolio/lgbouma_professional/static/movies"
cp $r/movie_EW_results_showvel_20260818/movie_EW_showvel_20231111_20231112.mov $blogdir/.
cd $blogdir
./convert_mov2mp4.sh
rm movie_EW_showvel_20231111_20231112.mov
mv $blogdir/movie_EW_showvel_20231111_20231112.mp4 $blogdir/LP_12-502_dbsp_Nov2023.mp4
cd $paperdir
cp $blogdir/LP_12-502_dbsp_Nov2023.mp4 figure_additions/.

# f6: FourStar vs TESS
cp $r/fourstar_vs_tess/fourstar_vs_tess_phase.pdf f6.pdf

# f7: MuSCAT vs TESS
cp $r/multicolor_phase/TIC262400835_multicolor_phase_stacked_reextracted.pdf f7.pdf

# f8: KeplerCam & Tierras vs TESS
cp $r/multicolor_phase/LP12502_multicolor_phase_stacked_Dec08Dec15.pdf f8.pdf

# f9: empirical power law fits to dip depth data
cp $r/beta_fit/combined/with_sigma_int/money_plot_beta_with_sigma_int.pdf f9.pdf

# f10: diluted transmission spectrum
cp $r/diluted_transmission_spectrum/diluted_transmission_spectrum.pdf f10.pdf

# f11: cartoon origin scenarios
cp $r/cgcd_cartoons/cgcd_cartoons_b_mixed_geom.pdf f11.pdf 

###############################################

# supp

# f13, f14: DBSP no dip nights
cp $r/EW_results/DBSP_HAlpha_mosaic_supp.pdf f13.pdf
cp $r/EW_results/LP12-502_Nov2024_dbsp_vs_tess.pdf f14.pdf

# f12: HIRES no dip nights
cp $r/hires_lp12502_mosaic/hires_lp12502_mosaic.pdf f12.pdf

# f15, f16: Dip fitting mosaics
cp $r/dip_mosaic/dip_fit_mosaic_1.pdf f15.pdf
cp $r/dip_mosaic/dip_fit_mosaic_2.pdf f16.pdf
