Title: Transient Corotating Clumps Around Adolescent Low-Mass Stars From Four Years of TESS

Authors: Luke G. Bouma, Rahul Jayaraman, Saul Rappaport, Luisa M. Rebull, Lynne A. Hillenbrand, Joshua N. Winn, Alexandre David-Uraz, and Gaspar A Bakos.

Contents: Column definitions for table1_MRT.csv

Defintiions are as follows.

ticid: TIC8.2 source identifer
dr2_source_id: Gaia DR2 source identifier (determined from TIC8.2).
dr3_source_id: Gaia DR3 source identifier (determined from dr2_source_id and gaiaedr3.dr2_neighbourhood table at Gaia archive)
angular_distance: Gaia DR2 -> DR3 crossmatch angular distance [arcseconds]
magnitude_difference: Gaia DR2 -> DR3 crossmatch G-band magnitude difference [mag]
dr3_ra: Gaia DR3 right ascension [deg]
dr3_dec: Gaia DR3 declination [deg]
tic8_Tmag: TIC8.2 TESS-band magnitude [mag]
dr3_dist_pc: Geometric Gaia DR3 distance [pc]
dr3_bp_rp: Gaia DR3 G_BP - G_RP color [mag]
dr3_ruwe: Gaia DR3 renormalized unit weight error
period: Adopted source period [hours]
assoc: Adopted source association
age: Adopted age corresponding to source association [Myr]
teff_sedfit: Effective temperature determined from SED fitting analysis [K]
teff_sedfit_perr: +1sigma statistical uncertainty on teff_sedfit [K]
teff_sedfit_merr: -1sigma statistical uncertainty on teff_sedfit [K]
rstar_sedfit: Stellar radius determined from SED fitting analysis [R_sun]
rstar_sedfit_perr: +1sigma statistical uncertainty on rstar_sedfit [R_sun]
rstar_sedfit_merr: -1sigma statistical uncertainty on rstar_sedfit [R_sun]
mass_parsec: Stellar mass determined from PARSECv1.2 model interpolation [M_sun]
mass_parsec_perr: +1sigma statistical uncertainty on mass_parsec [M_sun]
mass_parsec_merr: -1sigma statistical uncertainty on mass_parsec [M_sun]
Rcr_over_Rstar: (GM/Omega^2)^{1/3}, divided by rstar_sedfit
P2_hr: Secondary period present in TESS light curve, if present [hours]
quality: 1=CPV / 0=Candidate CPV / -1=Impostor
binarityflag:  Three-bit binarity flag, defined in end-notes of Table 1 in manuscript
N_sectors: number of of TESS sectors for which any data are expected to be acquired between 2018 July and 2024 Oct
