import numpy as np, pandas as pd, matplotlib.pyplot as plt
from os.path import join
import os
from complexrotators.paths import LOCALDIR, DATADIR, TABLEDIR
from aesthetic.plot import set_style, savefig
from numpy import array as nparr

def get_twomin():

    csvpath = "/Users/luke/local/SPOCLC/gaia_X_spoc2min_merge.csv"
    tdf = pd.read_csv(csvpath)

    sel = (
        (tdf.bp_rp > 1.5) &
        (tdf.TESSMAG < 16) &
        (tdf.parallax > 1e3*(1/150)) &
        (tdf.M_G > 4) &
        (tdf.SECTOR <= 55)
    )
    # S1-S55, sel fn applied
    stdf = tdf[sel]

    return stdf


def main():

    csvpath = 'twomin_banyan.csv'
    if not os.path.exists(csvpath):
        from cdips.utils.gaiaqueries import given_source_ids_get_gaia_data

        tdf = get_twomin()
        stdf = tdf.drop_duplicates('dr2_source_id', keep='first')

        source_ids = np.array(stdf.dr2_source_id).astype(np.int64)

        groupname = 'twomin_cpv_targetsample'
        gdf = given_source_ids_get_gaia_data(
            source_ids, groupname, n_max=65999, overwrite=False,
            enforce_all_sourceids_viable=True, savstr='',
            which_columns=(
                'g.source_id, g.ra, g.dec, g.parallax, g.pmra, g.pmdec, g.phot_g_mean_mag, '
                'g.pmra_error, g.pmdec_error, g.parallax_error'
            ),
            table_name='gaia_source', gaia_datarelease='gaiadr2',
            getdr2ruwe=False
        )
        gdf.to_csv(csvpath, index=False)

    else:
        gdf = pd.read_csv(csvpath)

    import sys
    sys.path.append("/Users/luke/Dropbox/proj/banyan_sigma")
    from core import membership_probability

    fn = lambda x: np.array(x).astype(float)

    ra, dec = fn(gdf.ra), fn(gdf.dec)

    pmra, pmdec = fn(gdf.pmra), fn(gdf.pmdec)

    epmra, epmdec = fn(gdf.pmra_error), fn(gdf.pmdec_error)

    plx, eplx = fn(gdf.parallax), fn(gdf.parallax_error)

    output = membership_probability(ra=ra, dec=dec, pmra=pmra, pmdec=pmdec,
                                    epmra=epmra, epmdec=epmdec, plx=plx, eplx=eplx,
                                    use_plx=True, use_rv=False)

    N_input = len(gdf)

    N_field_gt_95 = (output['ALL']['FIELD'] > 0.95).sum()
    N_field_gt_99 = (output['ALL']['FIELD'] > 0.99).sum()

    print(f'N_input: {N_input}')
    print(f'N_field_gt_95: {N_field_gt_95}')
    print(f'N_field_gt_99: {N_field_gt_99}')


if __name__ == "__main__":
    main()
