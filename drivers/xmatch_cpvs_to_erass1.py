"""
crossmatch eRASS1 coronal emitters against the CPV list from Bouma+2024
"""

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from complexrotators.paths import LITDIR, PAPERDIR
from os.path import join
from astropy import constants as c, units as u

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from aesthetic.plot import set_style, savefig


edf = pd.read_csv(join(LITDIR, "Freund_2024_eRASS1_likely_coronal.csv"))
edf = edf[edf.PLX > 0]
distances_pc = 1 / (edf.PLX * 1e-3)
Lx = (
    np.array(edf['Fx']) * u.erg * (u.s)**(-1) * (u.cm)**(-2) *
    4*np.pi*( np.array(distances_pc) * u.pc)**2
).to(u.erg/u.s).value
edf['Lx'] = Lx


b24_df = pd.read_csv(join(PAPERDIR, "table1_MRT.csv"), sep="|")
b24_df = b24_df[b24_df.quality >= 0]

mdf = b24_df.merge(edf, left_on="dr3_source_id", right_on="CTP_ID")

def get_f_G(mag):
    f_0 = 1
    mag_0 = 25.6884
    return f_0 * 10**( -0.4 * (mag - mag_0) )

lambda0 = 640.50 * u.nm
E = (c.h * c.c / lambda0).to(u.erg).value
A = ((0.7278*u.m**2).to(u.cm**2)).value
# estimate flux in erg/s/cm^2
f_G = E * get_f_G(mdf.G) / A
log10_Fx_by_Fg = np.log10(mdf.Fx/f_G)
mdf['log10_Fx_by_Fg'] = log10_Fx_by_Fg

# make scatter plot HRD

fig, ax = plt.subplots(figsize=(8,6))
set_style("clean")

sedf = edf[edf.PLX > 10]
ax.scatter(edf['BP_RP'], edf['Lx'], c='darkgray', s=1, zorder=1,
           linewidths=0, label='Freund+24 coronal, <100pc')

smdf = mdf[mdf.quality==1]
ax.scatter(smdf['BP_RP'], smdf['Lx'], c='k', s=8, label='B+24 CPV', zorder=5)
smdf = mdf[mdf.quality==0]
ax.scatter(smdf['BP_RP'], smdf['Lx'], c='darkgray', s=8, label='Maybe CPV',
           zorder=4)
smdf = mdf[mdf.binarityflag!=0]
ax.scatter(smdf['BP_RP'], smdf['Lx'], c='red', s=12, label='binary?', zorder=3)

ax.set_yscale('log')
ax.set_xlim([0,4.5])
ax.set_ylim([1e27, 1e32])
ax.legend(loc='best', fontsize='x-small')

ax.set_xlabel("BP-RP")
ax.set_ylabel("Lx")

outpath = 'eRASS1_CPV_HRD.png'
savefig(fig, outpath, writepdf=0)
