"""
Wed Nov  1 11:53:45 2023

table1_MRT.csv: updated using uniform SED fitting priors

table1_MRT.csv.BACKUP: original submitted version.
"""
import numpy as np, pandas as pd

df0 = pd.read_csv("table1_MRT.csv", sep="|")
df1 = pd.read_csv("table1_MRT.csv.BACKUP", sep="|")

cols = df0.columns

changedcols = ['teff_sedfit', 'teff_sedfit_perr', 'teff_sedfit_merr',
               'rstar_sedfit', 'rstar_sedfit_perr', 'rstar_sedfit_merr',
               'mass_parsec', 'mass_parsec_perr', 'mass_parsec_merr',
               'Rcr_over_Rstar', 'assoc', 'P2_hr']

for c in cols:
    if c in changedcols:
        continue
    print(f"{c}...")
    np.testing.assert_allclose(df0[c], df1[c])

pd.options.display.max_rows = 5000
DOPRINT = 0
if DOPRINT:
    print(pd.concat([df0['assoc'], df1['assoc']], axis=1))
    print(pd.concat([df0['P2_hr'], df1['P2_hr']], axis=1))

printcols = ['teff_sedfit', 'rstar_sedfit', 'mass_parsec', 'Rcr_over_Rstar',
             'teff_sedfit_perr', 'rstar_sedfit_perr', 'mass_parsec_perr']
for p in printcols:
    _df = pd.concat([df0[p], df1[p]], axis=1)
    _df['diff'] = df0[p] - df1[p]
    print(_df)
