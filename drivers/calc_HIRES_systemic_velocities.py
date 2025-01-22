import os
from glob import glob
from os.path import join
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from cdips_followup.spectools import get_naive_rv
from complexrotators.paths import DATADIR, RESULTSDIR

outdir = join(RESULTSDIR, 'HIRES_systemic_velocities')
if not os.path.exists(outdir): os.mkdir(outdir)

def main(RUNDICT, fitsdir, VBROAD, starid, run_in_parallel):

    rv_list = []
    drvs = []
    std_rvs = []
    fitsnames = []

    for starname,v in RUNDICT.items():

        fitsname, teff, logg, rv_expected = v
        fitsnames.append(fitsname)

        dirstarname = starname
        if "_" in starname:
            dirstarname = starname.split("_")[0]

        spectrum_path = join(fitsdir, fitsname)

        localdir = join(os.path.expanduser('~'), 'local')
        synth_path = join(
            localdir,
            f'synthetic_spectra/PHOENIX_MedRes/'
            f'lte0{teff:d}-{logg:.2f}-0.0'
            f'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
        )

        df = get_naive_rv(spectrum_path, synth_path, outdir, make_plot=1,
                          run_in_parallel=run_in_parallel, vbroad=VBROAD)

        # recalculate!
        rvs = np.array(df['rv_chisq_kms'])
        bc = np.array(df['bc_kms'])
        ZP = 81.527 # HIRES instrumental velocity ZP on r-chip is 81.527 +/- 0.664 km/s
        rvs = rvs - bc + ZP  # this is the systemic RV!!! (at the order level)

        # good for TIC402980664
        rchip_good_orders = [0,1,2,3,4,7,8,13,14,15] # default for ZP calc: [0,2,3,4,5,6,11]

        # TIC141146667
        rchip_good_orders = [2,6,13] # default for ZP calc: [0,2,3,4,5,6,11]

        sel_rvs = rvs[np.array(rchip_good_orders)]
        df['meangoodorder_rv_chisq_minus_bc_kms'] = np.round(np.nanmean(sel_rvs), 4)

        rv_std = np.round(np.nanstd(sel_rvs), 4)
        df['stdgoodorder_rv_chisq_minus_bc_kms'] = rv_std

        rv = df['meangoodorder_rv_chisq_minus_bc_kms'].iloc[0]

        if rv_expected is not None:
            drv = rv_expected - rv
            print(42*'*')
            print(f"{starname}: drv = {drv:.1f}")
            print(42*'*')
        else:
            drv = None

        rv_list.append(rv)
        std_rvs.append(rv_std)

        drvs.append(drv)

    df = pd.DataFrame({
        'fitsname': fitsnames,
        'starname': list(RUNDICT.keys()),
        'systemic_rv': rv_list,
        'std_rv_stat': std_rvs,
        'std_rv_zp': np.ones(len(rv_list))*0.664
    })
    df['std_systemic_rv'] = np.sqrt(
        df.std_rv_stat ** 2
        +
        df.std_rv_zp ** 2
    )

    outpath = join(outdir, f'{starid}_systemic_velocities.csv')
    df.to_csv(outpath, index=False)
    print(f'Wrote {outpath}')

    import IPython; IPython.embed()


if __name__ == "__main__":

    run_in_parallel = 1

    starid = 'TIC141146667'
    VBROAD = 130
    fitsdir = join(DATADIR, 'spectra/HIRES/TIC141146667_DEBLAZED')
    chip = 'r'
    fitspaths = np.sort(glob(join(fitsdir, f"{chip}j*fits")))
    RUNDICT = {}
    for fitspath in fitspaths:
        ind = os.path.basename(fitspath).rstrip(".fits").lstrip(f"{chip}")
        key = f'TIC141146667_{ind}'
        # file name, grid teff, %2f grid logg, expected RV (if there is one)
        val = [os.path.basename(fitspath), 3000, 4.00, 0]
        RUNDICT[key] = val


    ## file name, grid teff, %2f grid logg, expected RV (if there is one)
    #RUNDICT = {'TIC402980664': ['rj520.93.fits', 3300, 5.00, -12.33] }

    main(RUNDICT, fitsdir, VBROAD, starid, run_in_parallel)
