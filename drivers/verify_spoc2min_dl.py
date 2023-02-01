"""
Once the tesscurl_sector_XX_lc.sh scripts have been executed, you need to
verify that all available 2-minute light curves were actually downloaded
correctly by the CURL scripts.  In cases for which it failed, you need to
re-download them.  The verification step requires opening each of the files,
which is a bit slow (takes 4 minutes to open 20k light curves).

Contents:
    check_if_files_open
    make_curl_rectification_scripts
"""
from glob import glob
import os
from os.path import join
import numpy as np, pandas as pd
from astrobase.imageutils import get_header_keyword_list
from datetime import datetime
from complexrotators.paths import SPOCDIR

def check_if_files_open():

    sectors = range(1, 59)

    csvpaths = []

    for sector in sectors:

        csvpath = join(SPOCDIR, "FILE_OPENS", f"sector-{sector}_opens.csv")
        if os.path.exists(csvpath):
            print(f"Found {csvpath}")
            csvpaths.append(csvpath)
            continue

        fitspaths = glob(join(SPOCDIR, f"sector-{sector}", "*fits"))

        file_opens = []

        N = len(fitspaths)

        for ix, f in enumerate(fitspaths):

            if ix % 1000 == 0:
                print(f"{datetime.utcnow().isoformat()}: sector {sector}: {ix}/{N}...")

            it_opens = 1

            try:
                _ = get_header_keyword_list(f, ["TICID"], ext=0)
            except Exception as e:
                it_opens = 0

            file_opens.append(it_opens)

        df = pd.DataFrame({
            'fitspath': fitspaths,
            'file_opens': file_opens
        })
        df.to_csv(csvpath, index=False)
        print(f"Wrote {csvpath}")
        csvpaths.append(csvpath)

    return csvpaths

def make_curl_rectification_scripts(sectors):

    for sector in sectors:

        csvpath = join(SPOCDIR, "FILE_OPENS", f"sector-{sector}_opens.csv")
        assert os.path.exists(csvpath)

        df = pd.read_csv(csvpath)

        bad_sel = df.file_opens.astype(int) == 0

        bad_fitspaths = np.array(df[bad_sel].fitspath)

        for f in bad_fitspaths:
            os.remove(f)
            print(f"Deleted {f}")

        if len(bad_fitspaths) > 0:

            curl_lines = []
            for f in bad_fitspaths:
                _f = os.path.basename(f)
                l = (
                    f"curl -C - -L -o {_f} "
                    f"https://mast.stsci.edu/api/v0.1/Download/file/"
                    f"?uri=mast:TESS/product/{_f} \n"
                )
                curl_lines.append(l)

            outpath = join(
                SPOCDIR, "CURL_SCRIPTS", f"rectify_sector_{sector}.sh"
            )

            with open(outpath, 'w') as f:
                f.write('#!/bin/sh\n')
                f.writelines(curl_lines)

            print(f"Wrote {outpath}")


def main():

    sectors = range(1, 59)
    _ = check_if_files_open()
    make_curl_rectification_scripts(sectors)


if __name__ == "__main__":
    main()
