"""
Once verify_spoc2min_dl.py has been executed, this script just opens all of the
available *-s_lc.fits 2-minute cadence light curves, and pull all relevant FITS
header information into CSV files that are more readily accessible for
visualization and cross-matching purposes.
"""

import pandas as pd
from glob import glob
import os
from os.path import join
from astrobase.imageutils import get_header_keyword_list
from datetime import datetime

from complexrotators.paths import SPOCDIR

sectors = range(1, 59)

keyword_list = ["DATE", "TSTART", "TSTOP", "DATE-OBS", "DATE-END", "PROCVER",
                "OBJECT", "TICID", "SECTOR", "CAMERA", "CCD", "PXTABLE",
                "RA_OBJ", "DEC_OBJ", "PMRA", "PMDEC", "TESSMAG", "TEFF",
                "LOGG", "MH", "RADIUS", "TICVER", "CRMITEN", "CRBLKSZ",
                "CRSPOC", "CHECKSUM"]

for sector in sectors:

    csvpath = join(
        SPOCDIR, "METADATA", f"spoc2min_metadata_sector-{sector}.csv"
    )
    if os.path.exists(csvpath):
        print(f"Found {csvpath}")
        continue

    fitspaths = glob(join(SPOCDIR, f"sector-{sector}", "*fits"))

    rows = []
    N = len(fitspaths)

    for ix, fitspath in enumerate(fitspaths):

        if ix % 1000 == 0:
            print(f"{datetime.utcnow().isoformat()}: sector {sector}: {ix}/{N}...")

        try:
            d = get_header_keyword_list(fitspath, keyword_list, ext=0)
        except Exception as e:
            print(f"Got {e} for {fitspath}")
        d['fitspath'] = os.path.abspath(fitspath)
        rows.append(d)

    df = pd.DataFrame(rows)

    df.to_csv(csvpath, index=False)
    print(f"Wrote {csvpath}")
