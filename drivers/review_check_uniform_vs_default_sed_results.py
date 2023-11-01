import numpy as np, pandas as pd
import os
from os.path import join
from make_cpv_table import get_sedfit_results
from complexrotators.paths import PAPERDIR

csvpath = join(PAPERDIR, 'table1_MRT.csv')
df = pd.read_csv(csvpath, sep="|")

params = ['teff', 'rstar']

ticids = list(df.ticid)

dTs, sigTs, dRs, sigRs, sigTratio = [], [], [], [], []

ticids = np.array(ticids)


for ticid in ticids:

    uniform_sed_df = get_sedfit_results(ticid, uniformprior=1)

    dteff = int(
        df.loc[df.ticid == ticid, f'teff_sedfit'].iloc[0]
        -
        uniform_sed_df['teff_sedfit'].iloc[0]
    )
    sigteff = np.nanmean([
        df.loc[df.ticid == ticid, f'teff_sedfit_perr'],
        df.loc[df.ticid == ticid, f'teff_sedfit_merr'],
    ])
    u_sigteff = np.nanmean([
        uniform_sed_df['teff_sedfit_perr'].iloc[0],
        uniform_sed_df['teff_sedfit_merr'].iloc[0]
    ])

    sigTratio.append(
        u_sigteff / sigteff
    )
    dTs.append(dteff)
    sigTs.append(sigteff)

    dr = float(
        df.loc[df.ticid == ticid, f'rstar_sedfit'].iloc[0]
        -
        uniform_sed_df['rstar_sedfit'].iloc[0]
    )
    sigr = np.nanmean([
        df.loc[df.ticid == ticid, f'rstar_sedfit_perr'],
        df.loc[df.ticid == ticid, f'rstar_sedfit_merr'],
    ])
    dRs.append(dr)
    sigRs.append(sigr)

    print(f"{ticid}: ΔT={dteff}, σT={sigteff}")
    print(f"{ticid}: ΔRstar={dr:.3f}, σR={sigr:.3f}")


dTs = np.array(dTs)
sigTs = np.array(sigTs)
dRs = np.array(dRs)
sigRs = np.array(sigRs)
sigTratio = np.array(sigTratio)

print(f'N for which |T_u-T_orig| / σTorig > 1 = {(np.abs(dTs)/sigTs > 1).sum()}')
print(f'N for which |R_u-R_orig| / σRorig > 1 = {(np.abs(dRs)/sigRs > 1).sum()}')
print(f'N for which |T_u-T_orig| / σTorig > 2 = {(np.abs(dTs)/sigTs > 2).sum()}')
print(f'N for which |R_u-R_orig| / σRorig > 2 = {(np.abs(dRs)/sigRs > 2).sum()}')
print(f'N for which σT_u / σTorig > 1.3 = {(sigTratio > 1.5).sum()}')
print(f'N for which σT_u / σTorig < 0.7 = {(sigTratio < 0.8).sum()}')
print(f'Ntot = {len(dTs)}')

import IPython; IPython.embed()
#import matplotlib.pyplot as plt
#plt.hist(sigTratio, bins=30)
#plt.xlabel('σTeff uniform / σTeff "original"')
#plt.savefig('')
