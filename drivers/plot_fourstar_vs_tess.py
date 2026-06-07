import numpy as np, matplotlib.pyplot as plt, pandas as pd
import pickle
from os.path import join
from aesthetic.plot import set_style
import os

fourstardir = join(os.path.expanduser("~"),
                   "Dropbox/proj/cpv/data/photometry/FourStar_RAW_DATA/fourstar_20250210/PHOTOMETRY")
photpath = join(fourstardir, 'TIC300651846_aperture_photometry.xls')

t0_BJD_TDB = 2460693.243837336
period = 0.3439207894624283

fsdf = pd.read_csv(photpath, sep='\t')
# tcorr = 0.001415770035237074 # bjd corr (add this)
tcorr_minute = 60
tcorr_day = tcorr_minute / (60*24)
t = fsdf['BJD_TDB'] + tcorr_day # time: units days
f = fsdf['rel_flux_C1'] - 0.76
mask = (t - 2460717) > 0.87 # drop the last ~half hour
t, f = t[~mask], f[~mask]

csv_path = join(os.path.expanduser("~"),
                'Dropbox/proj/cpv/papers/Bouma_2026_cgcd/tables/phot_FourStar_20250210.csv')
pd.DataFrame({'t_BJD_TDB': t.values,
              'rel_flux': f.values - np.nanmean(f.values)}).to_csv(csv_path, index=False)

phase = (t - t0_BJD_TDB) / period - np.floor( (t-t0_BJD_TDB) / period )

# Time-bin FourStar flux at dt = 5 minutes
dt_days = 6/(24*60)
_t = np.array(t, dtype=float)
_f = np.array(f, dtype=float)
mask = np.isfinite(_t) & np.isfinite(_f)
_t, _f = _t[mask], _f[mask]
if len(_t):
    t0 = _t.min()
    bin_ix = np.floor((_t - t0)/dt_days).astype(int)
    df_fs = pd.DataFrame({
        'bin_ix': bin_ix,
        't': _t,
        'f': _f,
    })
    gb = df_fs.groupby('bin_ix')
    t_binned = gb['t'].mean().to_numpy()
    f_binned = gb['f'].mean().to_numpy()
    phase_binned = (t_binned - t0_BJD_TDB)/period
    phase_binned = phase_binned - np.floor(phase_binned)
else:
    t_binned = np.array([])
    f_binned = np.array([])
    phase_binned = np.array([])

tessdir = '/Users/luke/local/complexrotators/cpv_finding/spoc2min_debug'
tesspath = join(tessdir, '300651846_S0088_120sec_spoc2min_cpv_periodsearch.pkl')
with open(tesspath, 'rb') as __f:
    d = pickle.load(__f)
_t, _f = d['times'], d['fluxs']

_phase = (_t - t0_BJD_TDB) / period - np.floor( (_t - t0_BJD_TDB) / period )

# Phase-bin TESS flux to 200 points per cycle
nbins = 200
_t = np.array(_t, dtype=float)
_f = np.array(_f, dtype=float)
mask = np.isfinite(_t) & np.isfinite(_f)
_t, _f = _t[mask], _f[mask]
_phase = (_t - t0_BJD_TDB)/period
_phase = _phase - np.floor(_phase)
if len(_phase):
    p_ix = np.floor(_phase * nbins).astype(int)
    # Ensure last edge goes to last bin instead of index nbins
    p_ix[p_ix==nbins] = nbins-1
    df_tess = pd.DataFrame({
        'p_ix': p_ix,
        'phase': _phase,
        'flux': _f,
    })
    gbt = df_tess.groupby('p_ix')
    phase_bin_center = gbt['phase'].mean().to_numpy()
    flux_phase_binned = gbt['flux'].mean().to_numpy()
    # Sort by phase for plotting
    srt = np.argsort(phase_bin_center)
    phase_bin_center = phase_bin_center[srt]
    flux_phase_binned = flux_phase_binned[srt]
else:
    phase_bin_center = np.array([])
    flux_phase_binned = np.array([])


phase[phase > 0.5] -= 1
phase_binned[phase_binned > 0.5] -= 1
_phase[_phase > 0.5] -= 1
phase_bin_center[phase_bin_center > 0.5] -= 1

# Make plot overlaying raw and binned data
set_style("science")
plt.close('all')
fig, ax = plt.subplots(figsize=(2.,4))

# FourStar raw and time-binned, plotted vs phase
ax.scatter(phase, (f - 1)*100, s=1, color='maroon', alpha=0.15, marker='o', linewidths=0)
if len(phase_binned):
    # Sort for line plotting
    srt_fs = np.argsort(phase_binned)
    ax.scatter(phase_binned[srt_fs], (f_binned[srt_fs] - 1)*100, s=6, color='maroon',
               alpha=0.9, marker='o', linewidths=0)

# TESS raw and phase-binned
ax.scatter(_phase, (_f - 1)*100, s=0.08, color='C0', alpha=0.3, marker='o', linewidths=0)
if len(phase_bin_center):
    ax.scatter(phase_bin_center, (flux_phase_binned - 1)*100, s=5, color='C0',
               alpha=0.9, marker='o', linewidths=0)

ax.set_xlabel('Phase ($P$=8.26 hr)', fontsize=6)
ax.set_ylabel('Relative flux (%)', fontsize=6)
ax.tick_params(axis='both', which='major', labelsize=5.5)
ax.set_xlim(-0.55,0.55)
ax.set_ylim(-8, 11)
ax.set_yticks([-5, 0, 5, 10])
ax.set_title('TIC 300651846: 2.09$\mu$m', fontsize=7)

# Annotations for datasets
ax.text(0.49, 10, 'FourStar\n2025/02/10', color='maroon', ha='right',
        va='center', fontsize=4.5)
ax.text(0.49, -0.5, 'TESS\nSimult.', color='C0', ha='right', va='top',
        fontsize=4.5)

# Filter response inset
bandpass_dir = join(os.path.expanduser("~"),
                    'Dropbox/proj/cpv/data/photometry/svo_bandpasses')
nb209 = np.loadtxt(join(bandpass_dir, 'LCO_FourStar.NB209.dat'))
tess_bp = np.loadtxt(join(bandpass_dir, 'TESS_TESS.Red.dat'))

nb209_wave = nb209[:,0] * 1e-4   # Angstrom -> micron
nb209_trans = nb209[:,1] / nb209[:,1].max()
tess_wave = tess_bp[:,0] * 1e-4
tess_trans = tess_bp[:,1] / tess_bp[:,1].max()

axins = ax.inset_axes([0.025, 0.06, 0.38, 0.12])
axins.set_facecolor('white')
axins.fill_between(tess_wave, tess_trans, color='C0', alpha=0.8, lw=0)
axins.fill_between(nb209_wave, nb209_trans, color='maroon', alpha=0.9, lw=0)
axins.set_xlim(0.3, 2.3)
axins.set_ylim(0, 1.3)
axins.set_xticks([0.5, 1.0, 1.5, 2.0])
axins.set_xticklabels(['0.5', '1', '1.5', '2'], fontsize=4.5)
axins.minorticks_off()
axins.tick_params(axis='x', which='major', direction='in', pad=1, length=2,
                  bottom=True, top=False)
axins.tick_params(axis='y', which='both', left=False, right=False)
axins.set_yticks([])
for _spine in ['top', 'left', 'right']:
    axins.spines[_spine].set_visible(False)
axins.set_xlabel('μm', fontsize=5, labelpad=1)
axins.text(0.77, 1, 'TESS', color='C0', ha='center', va='bottom', fontsize=4)
axins.text(2.09, 1, 'FourStar\n2.09μm', color='maroon', ha='center',
           va='bottom', fontsize=4)

fig.tight_layout()
fig.savefig('../results/fourstar_vs_tess/fourstar_vs_tess_phase.png', dpi=300,
           bbox_inches='tight')
fig.savefig('../results/fourstar_vs_tess/fourstar_vs_tess_phase.pdf', dpi=300,
           bbox_inches='tight')
