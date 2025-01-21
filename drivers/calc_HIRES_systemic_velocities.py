import os
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from cdips_followup.spectools import get_naive_rv

outdir = '/Users/luke/Dropbox/proj/cpv/results/HIRES_systemic_velocities'
if not os.path.exists(outdir): os.mkdir(outdir)
run_in_parallel = 1

# NOTE: tested in /cdips_followup/tests/test_naive_rv_zeropoint.py

RUNDICT = {
    # file name, grid teff, %2f grid logg, expected RV (if there is one)
    'TIC402980664': ['rj520.93.fits', 3300, 5.00, -12.33], # 5 & 11 are off!!  13,14,15 all good
    'TIC141146667_537.166': ['rj537.166.fits', 3000, 4.50, 0],
    'TIC141146667_537.167': ['rj537.167.fits', 3000, 4.50, 0],
    'TIC141146667_500.288': ['rj500.288.fits', 3000, 4.50, 0],
}

def main():

    rv_list = []
    drvs = []
    std_rvs = []

    for starname,v in RUNDICT.items():

        fitsname, teff, logg, rv_expected = v

        dirstarname = starname
        if "_" in starname:
            dirstarname = starname.split("_")[0]

        spectrum_path = (f'/Users/luke/Dropbox/proj/cdips_followup/data/'
                         f'spectra/HIRES/{dirstarname}/{fitsname}')
        synth_path = (f'/Users/luke/local/synthetic_spectra/PHOENIX_MedRes/'
                      f'lte0{teff:d}-{logg:.2f}-0.0'
                      f'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

        df = get_naive_rv(spectrum_path, synth_path, outdir, make_plot=1,
                          run_in_parallel=run_in_parallel)

        # recalculate!
        rvs = np.array(df['rv_chisq_kms'])
        bc = np.array(df['bc_kms'])
        ZP = 81.527 # HIRES instrumental velocity ZP on r-chip is 81.527 +/- 0.664 km/s
        rvs = rvs - bc + ZP  # this is the systemic RV!!! (at the order level)

        # good for TIC402980664
        rchip_good_orders = [0,1,2,3,4,7,8,13,14,15] # default for ZP calc: [0,2,3,4,5,6,11]

        # TIC141146667
        rchip_good_orders = [2,6] # default for ZP calc: [0,2,3,4,5,6,11]

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
    import IPython; IPython.embed()

    #drvs = np.array(drvs)
    #std_rvs = np.array(std_rvs)

    #weights = 1.0 / (std_rvs**2)
    ## calculate zero-point, and its uncertainty, using a weighted sum based on
    ## each uncertainty (measured via order to order scatter)
    #wmean = np.sum(drvs * weights) / np.sum(weights)
    #wstd = np.sqrt(1.0 / np.sum(weights))

    #zeropoint = 1. * wmean  + ZP
    #zerpoint_unc = 1. * wstd

    #plt.close("all")
    #fig, ax = plt.subplots()
    #ax.errorbar(np.arange(len(drvs)), drvs, yerr=std_rvs, ls=":", c='k', markersize=2)
    #ax.set_xlabel("star index (G8 left, M4.5 right)")
    #ax.set_ylabel("expected RV - my RV - zeropoint")
    #title = f'zeropoint = {zeropoint:.3f} +/- {zerpoint_unc:.3f} km/s'
    #ax.set_title(title)
    #outpath = os.path.join(outdir, "zeropoint.png")
    #fig.savefig(outpath, bbox_inches='tight', dpi=300)

    #print(title)
    #print(f"saved {outpath}")


if __name__ == "__main__":
    main()
