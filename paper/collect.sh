resultsdir="/Users/luke/Dropbox/proj/cpv/results"
artdir="/Users/luke/Dropbox/proj/cpv/art"
addir="/Users/luke/Dropbox/Documents/proposals/2023_05_ADAP/latex/figs/"

# 3-panel mosaic and cartoons
cp $resultsdir/lc_mosaic/lc_mosaic_fav3.pdf f1a.pdf
cp $artdir/twomodels.pdf f1b.pdf

# completeness
cp $resultsdir/fraction_coverage/fraction_coverage.pdf f2.pdf

# fig3 is annotated_lp12502.key exported

# full LC mosaic (50 panel)
cp $resultsdir/full_lcmosaic/full_lcmosaic_showtitles_binarytitle_Ngoodsectors_tic8_Tmag.pdf f4.pdf

# catalog properties
cp $resultsdir/catalogscatter/catalogscatter_showmaybe_emphruwe.pdf f5.pdf

# LC evolution
cp $resultsdir/beforeafter_mosaic/beforeafter_mosaic_showtitles_tic8_Tmag.pdf f6.pdf

# TIC4029 segments
cp $resultsdir/tic4029_segments/tic4029_segments.pdf f7.pdf

# TIC4029 phase timegroups one orbit per panel
cp $resultsdir/phase_timegroups_mosaic/TIC_402980664/TIC_402980664_P18.5611_2min_phase_timegroups_mosaic_ymin-4.8_ymax3.pdf f8.pdf

# TIC4029 river
cp $resultsdir/river/TIC_402980664/TIC402980664_river_gist_stern_truncated_-1_64_manual_20230617_mask_v0_nterms2_vmin-4.00000_vmax1.00000.pdf f9a.pdf
cp $resultsdir/river/TIC_402980664/TIC402980664_river_gist_stern_truncated_248_315_manual_20230617_mask_v0_nterms2_vmin-4.00000_vmax1.00000.pdf f9b.pdf
cp $resultsdir/river/TIC_402980664/TIC402980664_river_gist_stern_truncated_1233_1265_manual_20230617_mask_v0_nterms2_vmin-4.00000_vmax1.00000.pdf f9c.pdf
cp $resultsdir/river/TIC_402980664/TIC402980664_river_gist_stern_truncated_1411_1481_manual_20230617_mask_v0_nterms2_vmin-4.00000_vmax1.00000.pdf f9d.pdf

# magnetic b star
cp $resultsdir/magnetic_bstar_comparison/magnetic_bstar_comparison_showtitles_simpleB.pdf f10a.pdf
cp $resultsdir/magnetic_bstar_comparison/magnetic_bstar_comparison_showtitles_complexB.pdf f10b.pdf
cp $resultsdir/magnetic_bstar_comparison/magnetic_bstar_comparison_showtitles_complexM.pdf f10c.pdf

# appendix
# TIC3006 phase timegroup
cp $resultsdir/phase_timegroups/TIC_300651846/cycle-1_to_630TIC_300651846_P8.2540_2min_phase_timegroups.pdf f11a.pdf
cp $resultsdir/phase_timegroups/TIC_300651846/cycle2290_to_2680TIC_300651846_P8.2540_2min_phase_timegroups.pdf f11b.pdf

# TIC3006 river
cp $resultsdir/river/tic_300651846/TIC300651846_river_seismic_-1_630_vmin-6.00000_vmax6.00000.pdf f12a.pdf
cp $resultsdir/river/tic_300651846/TIC300651846_river_seismic_2290_2680_vmin-8.00000_vmax8.00000.pdf f12b.pdf

# brightness/dist comp
cp $resultsdir/literaturecomp/gmag_vs_distance_CPVs_showlit.pdf f13.pdf

# specturm windows
cp $resultsdir/spectrum_windows/*402980664*pdf f14a.pdf
cp $resultsdir/spectrum_windows/*146539195*pdf f14b.pdf
cp $resultsdir/spectrum_windows/*264599508*pdf f14c.pdf
