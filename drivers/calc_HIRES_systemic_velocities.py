import os
from numpy import array as nparr
from astropy.io import fits
from glob import glob
from os.path import join
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from cdips_followup.spectools import get_naive_rv
from complexrotators.paths import DATADIR, RESULTSDIR

outdir = join(RESULTSDIR, 'HIRES_systemic_velocities')
if not os.path.exists(outdir): os.mkdir(outdir)

def main(RUNDICT, fitsdir, VBROAD, starid, chip, run_in_parallel,
         chip_good_orders=None):

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

        df = get_naive_rv(spectrum_path, synth_path, outdir, chip, make_plot=1,
                          run_in_parallel=run_in_parallel, vbroad=VBROAD)

        # recalculate!
        rvs = np.array(df['rv_chisq_kms'])
        bc = np.array(df['bc_kms'])
        # as determined in test_naive_rv_zeropoint.py
        ZEROPOINTS = {
            'r': 80.956, # HIRES instrumental velocity ZP on r-chip (+/- 0.664 km/s)
            'i': 80.956+5.443, # shift based on M dwarf TiO bands; makes scales consistent
        }
        ZP = ZEROPOINTS[chip]

        rvs = rvs - bc + ZP  # this is the systemic RV!!! (at the order level)

        # NOTE: fine-tuned to TIC141146667
        if chip == 'r' and chip_good_orders is None:
            # 2,6 tic1411 (+7?)
            chip_good_orders = [2,7] # default for ZP calc: [0,2,3,4,5,6,11]
        elif chip == 'i' and chip_good_orders is None:
            chip_good_orders = [8]

        sel_rvs = rvs[np.array(chip_good_orders)]
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
    df['std_systemic_rv'] = np.round(np.sqrt(
        df.std_rv_stat ** 2
        +
        df.std_rv_zp ** 2
    ), 4)

    sel = (np.abs(df.std_rv_stat) < 20) & (np.abs(df.systemic_rv) < 50)
    sdf = df[sel]
    weights = 1.0 / (sdf.std_systemic_rv**2)

    wmean = np.sum(sdf.systemic_rv * weights) / np.sum(weights)
    wstd = np.sqrt(1.0 / np.sum(weights))
    print(f'{starid} {chip} {chip_good_orders}: RV_sys = {wmean:.2f} +/- {wstd:.2f} km/s')

    cstr = "_".join([str(c) for c in chip_good_orders])
    outpath = join(outdir, f'{starid}_systemic_velocities_{chip}chip_orders_{cstr}.csv')
    df.to_csv(outpath, index=False)
    print(f'Wrote {outpath}')


if __name__ == "__main__":

    run_in_parallel = 1

    starid = 'TIC141146667'
    VBROAD = 130
    fitsdir = join(DATADIR, 'spectra/HIRES/TIC141146667_DEBLAZED')

    for chip, chip_good_orders in zip(
        ['r', 'i', 'r', 'r', 'r'],
        [[2,7], [8], [2], [7], [6]]
    ):

        fitspaths = np.sort(glob(join(fitsdir, f"{chip}j*fits")))

        RUNDICT = {}
        mjds = []
        for fitspath in fitspaths:

            hdul = fits.open(fitspath)
            hdr = hdul[0].header
            mjd = hdr['MJD']
            mjds.append(mjd)
            hdul.close()

            ind = os.path.basename(fitspath).rstrip(".fits").lstrip(f"{chip}")
            key = f'TIC141146667_{ind}'
            # file name, grid teff, %2f grid logg, expected RV (if there is one)
            val = [os.path.basename(fitspath), 3000, 4.50, 0]
            RUNDICT[key] = val

        main(RUNDICT, fitsdir, VBROAD, starid, chip, run_in_parallel,
             chip_good_orders=chip_good_orders)

    print("*"*42)
    # join to make relative RV table
    dfs = []

    plt.close("all")
    fig, ax = plt.subplots(figsize=(4,3))

    for ix, (chip, order) in enumerate(zip(
        ['i', 'r', 'r', 'r'],
        [[8], [2], [7], [6]]
    )):
        cstr = "_".join([str(c) for c in order])
        csvpath = join(outdir, f'{starid}_systemic_velocities_{chip}chip_orders_{cstr}.csv')
        _df = pd.read_csv(csvpath)

        # hack for TIC1411 to eliminate "bad" epochs, which seem mostly to be
        # caused by cosmic rays
        sel = (np.abs(_df.std_rv_stat) < 20) & (np.abs(_df.systemic_rv) < 50)

        _df.loc[~sel, 'systemic_rv'] = np.nan

        _df['rel_rv'] = _df['systemic_rv'] - np.nanmean(_df.systemic_rv)

        _df['mjd'] = mjds

        dfs.append(_df)

        ax.scatter(_df['mjd'].astype(float), _df['rel_rv'], c=f'C{ix}', s=5,
                    label=f'{chip}, order = {order}', alpha=0.8, linewidths=0)

    mdf = pd.concat(dfs)
    mjds = nparr(mjds).astype(float)
    mean_rel_rv = nparr(
        mdf.groupby('mjd')['rel_rv'].transform('mean').iloc[:len(_df)]
    )

    count = nparr(
        mdf.groupby('mjd')['rel_rv'].transform('count').iloc[:len(_df)]
    )

    std_rel_rv = nparr((
        mdf.groupby('mjd')['rel_rv'].transform('std')
        /
        np.sqrt(mdf.groupby('mjd')['rel_rv'].transform('count'))
    ).iloc[:len(_df)])

    weights = 1 / std_rel_rv**2
    wmean = np.sum(mean_rel_rv * weights) / np.sum(weights)

    outdf = pd.DataFrame({
        'mjd': mjds,
        'jd': 2400000.5 + mjds,
        'mean_rel_rv': np.round(mean_rel_rv - wmean, 2),
        'count': count,
        'std_rel_rv': np.round(std_rel_rv, 2),
    })

    ax.errorbar(outdf['mjd'], outdf['mean_rel_rv'],
                yerr=outdf['std_rel_rv'], c='k', marker='o', elinewidth=1,
                capsize=1, lw=0, mew=0.5, markersize=0)

    ax.legend(loc='best', fontsize='xx-small')
    ax.set_xlabel('time')
    ax.set_ylabel('rv km/s')

    outdir = join(RESULTSDIR, 'HIRES_systemic_velocities')
    outpath = join(outdir, f'{starid}_rel_rv.png')
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close("all")

    csvpath = join(outdir, f'{starid}_rel_rv.csv')
    outdf.to_csv(
        csvpath, index=False
    )
    print(f'wrote {csvpath}')

    outdir = '/Users/luke/Dropbox/proj/hotcold/paper/tables'
    csvpath = join(outdir, f'tab_{starid}_rel_rv.csv')
    outdf[['mjd','mean_rel_rv','std_rel_rv']].to_csv(
        csvpath, index=False
    )
    print(f'wrote {csvpath}')

    rounddict = {
        'mean_rel_rv': 2,
        'std_rel_rv': 2,
    }
    formatters = {}
    for k,v in rounddict.items():
        formatters[k] = lambda x: np.round(x, v)

    texpath = join(outdir, f'tab_{starid}_rel_rv.tex')
    outdf[['jd','mean_rel_rv','std_rel_rv']].to_latex(
        texpath, index=False, formatters=formatters
    )
    print(f'wrote {texpath}')

    pd.options.display.max_rows=100
