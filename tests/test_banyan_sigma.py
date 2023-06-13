import numpy as np, pandas as pd
import sys

# implicitly, from banyan-sigma core import...
sys.path.append("/Users/luke/Dropbox/proj/banyan_sigma")
from core import membership_probability

ra=311.2911826481039
dec=-31.3425000799281

pmra=281.319
epmra=0.022
pmdec=-360.148
epmdec=0.019

plx=102.943
eplx=0.023

rv=-5.2
erv=0.7

output = membership_probability(ra=ra, dec=dec, pmra=pmra, pmdec=pmdec,
                                epmra=epmra, epmdec=epmdec, plx=plx, eplx=eplx,
                                rv=rv, erv=erv, use_plx=True, use_rv=True)

probs = np.array(output['ALL'].iloc[0].round(6))
assocs = np.array(output['ALL'].iloc[0].index)
banyan_df = pd.DataFrame({"assoc": assocs, "prob": probs})

print(banyan_df.sort_values(by='prob', ascending=False).head(5))

output = membership_probability(ra=ra, dec=dec, pmra=pmra, pmdec=pmdec,
                                epmra=epmra, epmdec=epmdec, plx=plx, eplx=eplx,
                                use_plx=True, use_rv=False)

probs = np.array(output['ALL'].iloc[0].round(6))
assocs = np.array(output['ALL'].iloc[0].index)
banyan_df = pd.DataFrame({"assoc": assocs, "prob": probs})

print(banyan_df.sort_values(by='prob', ascending=False).head(5))
