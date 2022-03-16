import os
import pandas as pd
from glob import glob
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR, TARGETSDIR

# begin options
# PLOTDIRNAME = 'Rahul_22b_phase'
# TARGETLISTNAME = '20220115_RahulJayararaman_CR_list.csv'
# PLOTDIRNAME = 'LAH_S43_to_S45_CRs'
# TARGETLISTNAME = 'LAH_TESS_GI_Sector_43_thru_45_CRs.csv'
PLOTDIRNAME = 'year3_LGB_ddt_phase'
TARGETLISTNAME = 'complex_rotators_ddt_merge_tic8.csv'

# end options

PLOTDIR = os.path.join(RESULTSDIR, PLOTDIRNAME)
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

targetlist = os.path.join(TARGETSDIR, TARGETLISTNAME)
df = pd.read_csv(targetlist)

manualperiodfile = os.path.join(TARGETSDIR, 'ticids_manual_periods.csv')
mpdf = pd.read_csv(manualperiodfile)

for ticid in df.ticid.astype(str):
    outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    row = mpdf.loc[mpdf.ticid.astype(str) == ticid]
    if len(row) > 0:
        manual_period = float(row.manual_period)
    else:
        manual_period = None

    phasepaths = glob(os.path.join(outdir, '*_phase.png'))

    if len(phasepaths) > 0:
        print(f"Found {phasepaths}, skipping.")
        continue

    print(f'Beginning {ticid}')

    try:
        cp.plot_phase(
            outdir,
            ticid='TIC_'+ticid,
            lc_cadences='2min_20sec',
            binsize_minutes=10,
            manual_period=manual_period
        )
    except Exception as e:
        print(f'ERROR! {ticid} failed.\n Reason was...{e}')
