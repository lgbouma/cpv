import os
from astropy.io import fits
from glob import glob
from os.path import join
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from complexrotators.paths import RESULTSDIR
from astropy.time import Time

def plot_ew_timeseries(
    linename = 'HGamma',
    linekey = 'Hγ',
    datestr = '20231112'
):
    fitsdir = f'/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX/{datestr}/p200_dbsp_blue_A/Science'
    ewdir = f'/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX/{datestr}/p200_dbsp_blue_A/Balmer_EWs'
    fitspaths = np.sort(glob(join(fitsdir, "spec1d_blue*LP_12-502*fits")))

    print(len(fitspaths))

    times, keys = [], []
    for f in fitspaths:
        hl = fits.open(f)
        time = hl[0].header['MJD']
        key = os.path.basename(f).replace(".fits","")
        times.append(time)
        keys.append(key)
        hl.close()

    csvpaths = [glob(join(ewdir, '*'+linename+"*"+key+'*_results.csv'))[0] for key in keys]

    dfs = [pd.read_csv(f) for f in csvpaths]

    fitted_ews = [df['Fitted_EW_mA'].iloc[0] for df in dfs] # milliangstr
    perr = [df['Fitted_EW_mA_perr'].iloc[0] for df in dfs]
    merr = [df['Fitted_EW_mA_merr'].iloc[0] for df in dfs]
    errs = np.array([merr,perr])

    times = np.array(times)+2400000.5 # mjd to jd

    t = Time(times, format='jd')
    # # then t.jd can be passed and
    # # manually ran thru https://astroutils.astronomy.osu.edu/time/utc2bjd.php
    # # -->learn it is a five minute difference(...not important, probably)
    # bjd_times = np.array([
    #     2460259.599272287, 2460259.608587248, 2460259.615961203, 2460259.623156718,
    #     2460259.630352222, 2460259.637547795, 2460259.644743333, 2460259.651938917,
    #     2460259.659134478, 2460259.666330027, 2460259.673525519, 2460259.680721057,
    #     2460259.688240211, 2460259.695435748, 2460259.702631262, 2460259.709826811,
    #     2460259.717022301, 2460259.724217827, 2460259.731413306, 2460259.738608912,
    #     2460259.745804483, 2460259.753000077, 2460259.760195556, 2460259.767391069,
    #     2460259.777480704, 2460259.784676286, 2460259.791871822, 2460259.799067427,
    #     2460259.806263033, 2460259.813458591, 2460259.820654150, 2460259.827849743,
    #     2460259.835045221, 2460259.842240733, 2460259.849436256, 2460259.856631756,
    #     2460259.863827280, 2460259.871022815, 2460259.878218303, 2460259.885413757,
    #     2460259.892609234, 2460259.899804768, 2460259.907040118, 2460259.914235664,
    #     2460259.921431198, 2460259.928626709, 2460259.935822221, 2460259.943017731,
    #     2460259.950213219, 2460259.957408811, 2460259.964604322, 2460259.971799878,
    #     2460259.978995435, 2460259.986191039, 2460259.993386503, 2460260.000582014
    # ])

    from aesthetic.plot import set_style, savefig
    set_style('clean')
    fig, ax = plt.subplots()
    ax.errorbar(times - 2460259, fitted_ews, yerr=0.75*errs)
    ax.set_ylabel(f'{linekey} EW [mA]')
    ax.set_xlabel('JD_UTC - 2460259')
    outdir = join(RESULTSDIR, 'DBSP_results')
    if not os.path.exists(outdir): os.mkdir(outdir)
    outpath = join(outdir, f'{linename}_vs_time_{datestr}.png')
    savefig(fig, outpath)

if __name__ == "__main__":

    plot_ew_timeseries(linename='HBeta', linekey='Hβ', datestr='20231111')
    plot_ew_timeseries(linename='HBeta', linekey='Hβ', datestr='20231112')

    plot_ew_timeseries(linename='HGamma', linekey='Hγ', datestr='20231111')
    plot_ew_timeseries(linename='HGamma', linekey='Hγ', datestr='20231112')
