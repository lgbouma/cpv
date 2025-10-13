import numpy as np, matplotlib.pyplot as plt, pandas as pd
import pickle
from os.path import join
from aesthetic.plot import set_style

fourstardir = '/Users/luke/Dropbox/proj/cpv/data/photometry/FourStar_RAW_DATA/fourstar_20250210/PHOTOMETRY'
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
fig, ax = plt.subplots(figsize=(2.2,4))

# FourStar raw and time-binned, plotted vs phase
ax.scatter(phase, f, s=1, color='maroon', alpha=0.15, marker='o', linewidths=0)
if len(phase_binned):
    # Sort for line plotting
    srt_fs = np.argsort(phase_binned)
    ax.scatter(phase_binned[srt_fs], f_binned[srt_fs], s=6, color='maroon',
               alpha=0.9, marker='o', linewidths=0)

# TESS raw and phase-binned
ax.scatter(_phase, _f, s=0.08, color='C0', alpha=0.3, marker='o', linewidths=0)
if len(phase_bin_center):
    ax.scatter(phase_bin_center, flux_phase_binned, s=5, color='C0',
               alpha=0.9, marker='o', linewidths=0)

ax.set_xlabel('Phase')
ax.set_ylabel('Relative flux')
ax.set_xlim(-0.55,0.55)
ax.set_ylim(0.94, 1.11)
ax.set_yticks([0.95, 1, 1.05, 1.1])
ax.set_title('TIC 300651846 (P=8.3 hr)', fontsize=7)

# Annotations for datasets
ax.text(-0.49, 1.10, 'FourStar 2.09μm narrow\n2025/02/10', color='maroon', ha='left',
        va='center', fontsize=5)
ax.text(-0.49, 0.96, 'TESS 0.6-1μm\nAvg 01/14-02/11', color='C0', ha='left', va='top',
        fontsize=5)

fig.tight_layout()
fig.savefig('../results/fourstar_vs_tess/fourstar_vs_tess_phase.png', dpi=300,
           bbox_inches='tight')
fig.savefig('../results/fourstar_vs_tess/fourstar_vs_tess_phase.pdf', dpi=300,
           bbox_inches='tight')
