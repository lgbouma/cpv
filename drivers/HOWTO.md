----------
A) For getting all the 2-minute SPOC data and its Gaia xmatch
----------

1. Get all the bulk download scripts locally (done on wh1).
2. Run them.
3. Verify their output using `verify_spoc2min_dl.py`
4. Transform the metadata to CSV using `get_spoc2min_metadata.py`
   This yields `MERGED_spoc2min_sector1_to_sector58.csv`.
5. Get the unique TICIDs from that with an in-line script.
6. Upload the resulting CSV file to CasJobs at MAST, but be sure to move the
   largest TICIDs to the top of the list so that it guesses the correct dtype.
7. Run a job like
  ```
  SELECT uticid, ID, GAIA into mydb.MyTable from MyDB.s1s58
  JOIN CatalogRecord ON uticid = CatalogRecord.ID
  ```
  within the `TESS_82` context.
8. Move output to `~/local/SPOCLC/MyTable_lukebouma_s1_s58_CasJobs_TIC8.csv`
9. Run `get_spoc2min_gaiainfo.py`
10. Run a job like
  ```
  select u.dr2_source_id, g.source_id, g.ra, g.dec, g.parallax,
  g.parallax_over_error, g.pmra, g.pmdec, g.phot_g_mean_mag, g.phot_bp_mean_mag,
  g.phot_rp_mean_mag, g.phot_g_mean_flux_over_error,
  g.phot_bp_mean_flux_over_error, g.phot_rp_mean_flux_over_error, g.bp_rp,
  g.g_rp, g.radial_velocity, g.l, g.b
  from user_lbouma.s1s58 as u, gaiadr2.gaia_source as g
  where u.dr2_source_id=g.source_id
  ```
  on the Gaia archive.
11. Rerun `get_spoc2min_gaiainfo.py`
12. Run `merge_spoc2min_metadata.py`
13. Run `get_gaia_X_spoc2min_merge.py

----------
B) To run the dip-based CPV-finder
----------

1. Go to find_CPVs.py.  Define your sample shell, e.g., "85to95pc_mkdwarf", and
  run it.  (NB: for the first paper, I did this on phtess3)

----------
C) To assess what TESS data exists for a set of TIC IDs
----------

1. Run assess_tess_holdings.py

----------
D) To produce the table for the catalog manuscript
----------
1. Run (A) then (B)
2. Manual label dip-based pipeline results w/ TagSpaces, at /results/cpvvetter/
3. Get list of TIC ID's from Rahul, move to /data/targetlists/
4. Manual label Rahul's CPV candidates, also at /results/cpvvetter/
5. Merge labels using `merge_lgb_rahul_classification_lists.py`
6. Run find_CPVs.py using sample_id=="2023catalog_LGB_RJ_concat", which
   generates all the vetting reports in a homogeneous way for the ms.
7. Run `make_cpv_table.py` (relies on logs from step #6).

----------
E) To run the SED analysis
----------
1. Download BT-Settl spectra from
  http://svo2.cab.inta-csic.es/theory/newov2/index.php
  BT-Settl - Theoretical Spectra
  -> I did Teff from 2600-5500K, logg 3.5-6, met 0
  -> BT-Settl/AGSS2009

2. Install astroARIADNE into a dedicate conda environment.  Hack
  astroARIADNE/plotter.py around line 1300 to read in the SVO subset of
  BT-Settl spectra in the ascii format.

  (The alternative path, of going to
  https://phoenix.ens-lyon.fr/Grids/BT-Settl/AGSS2009/SPECTRA/ ->
  BT-Settl_M-0.0_cool.tar and BT-Settl_M-0.0_a+0.0_hot.tar produces a separate
  weird format, which is also not what is wanted.  The FITS format for AGSS2009
  seems to not be available).

3. Test on one star using tests/test_ariadne_fitter.py (and
   test_ariadne_plots.py)

4. Run `run_SED_analysis.py`, in the py38_ariadne environment.

  * Inspect results as /results/ariadne_sed_fitting, to flag photometric
    outliers and stars with IR excesses.
  * Rerun any necessary cases, deleting the entire directory in each case.

5. Rerun `make_cpv_table.py` to incorporate results into table.

----------
F) To run the custom TIC4029 analysis
----------
1. Run `plot_tic4029_segments.py`, which writes times/fluxs/cadenceno's into
   /results/tic4029_segments.
2. Manually flag fine times in glue; write to the same directory.
3. Run `build_4029_mask.py`, with settings flagged to generate the mask you
   want, and the model of the desired complexity.
4. Make river plots (plot_river.py)
5. Make phased timegroup mosaic (plot_phase_timegroups_mosaic.py)

----------
G) To make the plots and values files for the manuscript
----------
1. Run `get_selection_function_values.py`

----------
H) To reduce WINERED spectra
----------
1. From ~/proj/WARP, run exec.sh in the py37_warp environment.
  - Tune individual frame id's, and obj id's, as required.

----------
I) For getting all the QLP data and their Gaia xmatch
----------

1. Get all the bulk download scripts locally (done on wh1).
2. Run them.
3. Use ~/QLP/get_QLP_ticid_lc_location.py to get CSV files of ticids and paths.
4. Get the unique TICIDs from each an in-line script.
5. Upload the resulting CSV file to CasJobs at MAST, but be sure to move the
   largest TICIDs to the top of the list so that it guesses the correct dtype.
6. Run a job like
  ```
  SELECT ticid, ID, GAIA into mydb.MyTable_v2 from MyDB.s1s55_qlp
  JOIN CatalogRecord ON ticid = CatalogRecord.ID
  ```
  within the `TESS_82` context.
7. Download output via MyDB -> MyTable_v2
8. Move output to `~/local/QLP/MyTable_v2_lukebouma_s1s55_qlp_CasJobs_TIC8.csv`
9. Run `get_QLP_gaiainfo.py`.  Upload the result to Gaia archive by _splitting
   it in half_, because the max number of Gaia DR2 source_id's in an upload is
   ~1e7.
10. Run a job like
  ```
  select u.dr2_source_id, g.source_id, g.ra, g.dec, g.parallax,
  g.parallax_over_error, g.pmra, g.pmdec, g.phot_g_mean_mag, g.phot_bp_mean_mag,
  g.phot_rp_mean_mag, g.phot_g_mean_flux_over_error,
  g.phot_bp_mean_flux_over_error, g.phot_rp_mean_flux_over_error, g.bp_rp,
  g.g_rp, g.radial_velocity, g.l, g.b
  from user_lbouma.qlp_s1s55_firsthalf as u, gaiadr2.gaia_source as g
  where u.dr2_source_id=g.source_id
  ```
  on the Gaia archive.
11. Move resulting .vot.gz files to ~/local/QLP
12. Rerun `get_QLP_gaiainfo.py`
13. Run `get_QLP_closesubset.py` to do some parallax cuts, to make more
    manageable.

----------
X) To search for transiting planets
----------
...todo
