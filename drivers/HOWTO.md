For getting all the 2-minute SPOC data and its Gaia xmatch
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
  from user_lbouma.s1s58 import as u, gaiadr2.gaia_source as g
  where u.dr2_source_id=g.source_id
  ```
  on the Gaia archive.
11. Rerun `get_spoc2min_gaiainfo.py`
12. Run `merge_spoc2min_metadata.py`
13. Run `get_gaia_X_spoc2min_merge.py

To run the CPV-finder
----------

1. Go to find_CPVs.py.  Define your sample shell, e.g., "85to95pc_mkdwarf", and
  run it.

To assess what TESS data exists for a set of TIC IDs
----------

1. Run assess_tess_holdings.py

To search for transiting planets
----------
(todo)
