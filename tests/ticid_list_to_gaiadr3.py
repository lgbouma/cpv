"""
Given a list of TICIDs, get the Gaia DR3 source_ids...
    and then see which have RVS spectra, and download the spectra.
"""
# options:
# tic -> gaia dr2 -> convert all to dr3

from astrobase.services.identifiers import tic_to_gaiadr2

from complexrotators.paths import TARGETSDIR
from os.path import join
import pandas as pd, numpy as np

csvpath = join(TARGETSDIR, "20230411_good_CPV_ticids_d_lt_150pc.csv")
df = pd.read_csv(csvpath)

dr2_source_ids = []
for ticid in np.array(df.ticid):
    dr2_source_ids.append(tic_to_gaiadr2(str(ticid)))

df['dr2_source_id'] = dr2_source_ids

from cdips.utils.gaiaqueries import given_dr2_sourceids_get_edr3_xmatch

gdf = given_dr2_sourceids_get_edr3_xmatch(
    np.array(dr2_source_ids).astype(np.int64), '20230411_good_CPVs',
    overwrite=0
)

# has DR2 <-> DR3 matches for this case
# NOTE: generally you will want to message around with the angular_distance and
# absolute_magnitude_difference keys to see what the correct choice is here.
sdf = gdf[gdf.abs_magnitude_difference < 1]

selcols = (
    'dr2_source_id dr3_source_id angular_distance magnitude_difference'.split()
)
sdf = sdf[selcols]

sdf['dr2_source_id'] = sdf.dr2_source_id.astype(str)
df['dr2_source_id'] = df['dr2_source_id'].astype(str)

mdf = df.merge(sdf, how='inner', on='dr2_source_id')
assert len(mdf) == len(df)

outpath = join(
    TARGETSDIR,
    "20230411_good_CPV_ticids_d_lt_150pc_X_DR2_X_DR3_sourceids.csv"
)
mdf.to_csv(outpath, index=False)
print(f'wrote {outpath}')
#####################################
# got the DR3 source_ids
#####################################

#####################################
# check which CPVs have RVS spectra #
#####################################
from cdips.utils.gaiaqueries import given_source_ids_get_gaia_data

groupname = '20230411_good_CPV_dlt150_dr3_lite'
dr3_source_ids = np.array(mdf.dr3_source_id).astype(np.int64)
gdf = given_source_ids_get_gaia_data(dr3_source_ids, groupname, n_max=10000,
                                     overwrite=False,
                                     enforce_all_sourceids_viable=True,
                                     savstr='', which_columns='*',
                                     table_name='gaia_source_lite',
                                     gaia_datarelease='gaiadr3')

sgdf = gdf[gdf['has_rvs']]
sgdf['dr3_source_id'] = sgdf['source_id'].astype(str)
mdf['dr3_source_id'] = mdf.dr3_source_id.astype(str)
smgdf = mdf.merge(sgdf, on='dr3_source_id', how='inner')

outpath = join(
    TARGETSDIR,
    "20230411_good_CPV_ticids_d_lt_150pc_X_DR2_X_DR3_allinfo_has_RVS.csv"
)
smgdf.to_csv(outpath, index=False)
print(f'wrote {outpath}')

from cdips.utils.gaiaqueries import given_dr3_sourceids_get_rvs_spectra
dr3_source_ids = np.array(smgdf['dr3_source_id']).astype(np.int64)
cache_id = '20230411_good_CPV_dlt150pc_hasRVS'
rvspaths = given_dr3_sourceids_get_rvs_spectra(dr3_source_ids, cache_id)
