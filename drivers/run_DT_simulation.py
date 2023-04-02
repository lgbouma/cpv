"""
Simulate doppler tomography recovery of CPV signals
"""
import numpy as np, matplotlib.pyplot as plt
from astropy import units as u, constants as c
from aesthetic.plot import set_style, savefig
from complexrotators.paths import RESULTSDIR
import os
from os.path import join

outdir = join(RESULTSDIR, "DT_simulation")
if not os.path.exists(outdir): os.mkdir(outdir)

##########################################
# BEGIN CONFIG PARAMS
star_id = "ID-2"
Vmag = 13.9
period = 7.9*u.hour
#t_exp = (10*u.min).to(u.day).value
t_exp = (12*u.min).to(u.day).value
t_read = (1*u.min).to(u.day).value
t_baseline = (7*u.hour).to(u.day).value  # at airmass < 2.5
A_spot = 0.035/2 # spot-induced semi-amplitude
dip_dict = {
    # (key) index: (values) phase_location, depth_in_percent, duration_phase_units
    'Dip0': {'loc': -0.35, 'depth':1.6, 'dur':0.18},
    'Dip1': {'loc': 0, 'depth':1.5,   'dur':0.198},
    'Dip2': {'loc': 0.33, 'depth':2, 'dur':0.195},
}
#ID-1
#dip_dict = {
#    # (key) index: (values) phase_location, depth_in_percent, duration_phase_units
#    'Dip0': {'loc': -0.35, 'depth':0.7, 'dur':0.198},
#    'Dip1': {'loc': -0.15, 'depth':2,   'dur':0.198},
#    'Dip2': {'loc': 0.025, 'depth':1.5,     'dur':0.15},
#}
# END CONFIG PARAMS
##########################################

# calculate number of exposures
N_exp = int(t_baseline / (t_exp + t_read))

# estimate vsini
Rstar = 1*u.Rsun
sini = 1
v = (2*np.pi*Rstar / period).to(u.km/u.s)
vsini = (v * sini).to(u.km/u.s).value

# estimate spot-induced signal
t0 = 0
N_theory = 1000
t = np.linspace(t0, t0+t_baseline, N_theory)
t_obs = np.linspace(t0, t0+t_baseline, N_exp)

P = period.to(u.day).value
phase_obs = (t_obs - t0)/P - np.floor((t_obs - t0)/P)
phase = (t - t0)/P - np.floor((t - t0)/P)

t0_spot = t0 + 0.*P

# inject spot signal at 2*Ω_rot
RV_spot = A_spot * vsini * np.sin( 2 * 2*np.pi*(t-t0_spot)/P )
RV_spot_obs = A_spot * vsini * np.sin( 2 * 2*np.pi*(t_obs-t0_spot)/P )

# estimate transient feature signals.  assume that in velocity they are
# basically "extra sinusoids", similar to the RM signal.
RV_transient = np.zeros(N_theory)
RV_transient_obs = np.zeros(N_exp)

for dip_ix, dip_vals in dip_dict.items():

    dip_phase, dip_depth_pct, dip_dur_phase = (
        dip_vals['loc'], dip_vals['depth'], dip_vals['dur']
    )

    t0_dip = 0.40*P + dip_phase*P

    dip_dur_time = dip_dur_phase*P

    A_dip = dip_depth_pct*1e-2

    # theoretical one...
    RV_dip = -1.* A_dip * vsini * np.sin( 2*np.pi*(t-t0_dip) / dip_dur_time )

    t_ing = t0_dip - dip_dur_time/2
    t_egr = t0_dip + dip_dur_time/2

    dip_dict[dip_ix]['t_ing'] = t_ing
    dip_dict[dip_ix]['t_egr'] = t_egr

    in_dip = ( t < t_egr ) & ( t > t_ing )
    RV_dip[~in_dip] = 0
    RV_transient += RV_dip

    # observed one... [would be better to do exposure time integration]
    RV_dip_obs = -1.* A_dip * vsini * np.sin( 2*np.pi*(t_obs-t0_dip) / dip_dur_time )
    in_dip_obs = ( t_obs < t_egr ) & ( t_obs > t_ing )
    RV_dip_obs[~in_dip_obs] = 0
    RV_transient_obs += RV_dip_obs

err = 0 # leave as zero, b/c we won't be measuring velocities anyway
RV_model = RV_spot + RV_transient
RV_obs = RV_spot_obs + RV_transient_obs + err

# make rv vs phase plot
set_style("clean")
fig, ax = plt.subplots()
ax.plot(
    phase, RV_model, c='k', lw=1
)
ax.scatter(
    phase_obs, RV_obs, c='k', marker='X', s=2
)
ax.update({
    'xlabel':'phase', 'ylabel':'RV [km/s]'
})
outpath = join(outdir, f"20230331_{star_id}_DT_simulation_rv_phase.png")
savefig(fig, outpath)

##########
# proceed to the line profile style analysis

# velocity resolution
R = 30000
dv = (c.c / R).to(u.km/u.s).value
velocity_grid = np.arange(-3*vsini, +3*vsini, dv)

from scipy.stats import norm
ccf_mean_model = []

sigma = vsini/(2**0.5)

factor = 1 / (sigma * np.sqrt(2*np.pi))

#
# N_velocity_elements X N_times
#
RV_zero = np.zeros_like(RV_model)
RV_zero_obs = np.zeros_like(RV_obs)

ccf_model = factor * np.exp(
    -0.5 * ((velocity_grid[:,None] - RV_zero[None,:])/sigma)**2
)
# NOTE: no noise yet in observations...
ccf_obs = factor * np.exp(
    -0.5 * ((velocity_grid[:,None] - RV_zero_obs[None,:])/sigma)**2
)

# get the sinusoidal spot model (in the residual)

# inject spot signal at 2*Ω_rot
def fn_ccf_resid_spot_model_fixed_time(A_spot, vsini, t0_spot, t, P, dv):
    RV = vsini * np.sin( 2 * 2*np.pi*(t-t0_spot)/P )
    sigma = 5*dv
    #factor = 1 / (sigma * np.sqrt(2*np.pi))
    factor = 1
    ccf_resid = - A_spot * factor * np.exp(
        -0.5 * ((velocity_grid - RV)/sigma)**2
    )
    return ccf_resid

def fn_ccf_resid_RM_model_fixed_time(A_dip, vsini, t0_dip, dip_dur_time, t, P, dv):

    t_ing = t0_dip - dip_dur_time/2
    t_egr = t0_dip + dip_dur_time/2

    in_dip = ( t < t_egr ) & ( t > t_ing )

    RV_dip = 1.* vsini * np.sin(0.3 * 2*np.pi*(t-t0_dip) / dip_dur_time )

    sigma = 2*dv

    #factor = 1 / (sigma * np.sqrt(2*np.pi))
    factor = 1
    ccf_resid = - A_dip * factor * np.exp(
        -0.5 * ((velocity_grid - RV_dip)/sigma)**2
    )

    if not in_dip:
        return 0*ccf_resid

    else:
        return ccf_resid

#
# spot ccfs: model and obs
#
ccf_resid_spot_model = []
for _t in t:
    ccf_resid_spot_model.append(
        fn_ccf_resid_spot_model_fixed_time(A_spot, vsini, t0_spot, _t, P, dv)
    )
ccf_resid_spot_model = np.array(ccf_resid_spot_model)

ccf_resid_spot_obs = []
for _t in t_obs:
    ccf_resid_spot_obs.append(
        fn_ccf_resid_spot_model_fixed_time(A_spot, vsini, t0_spot, _t, P, dv)
    )
ccf_resid_spot_obs = np.array(ccf_resid_spot_obs)

#
# dip ccfs: model and obs
#
ccf_resid_dips_model = np.zeros_like(ccf_resid_spot_model)
ccf_resid_dips_obs = np.zeros_like(ccf_resid_spot_obs)

for ix, (dip_ix, dip_vals) in enumerate(dip_dict.items()):

    ccf_resid_dip_model = []
    ccf_resid_dip_obs = []

    dip_phase, dip_depth_pct, dip_dur_phase = (
        dip_vals['loc'], dip_vals['depth'], dip_vals['dur']
    )
    A_dip = dip_depth_pct*1e-2
    t0_dip = 0.40*P + dip_phase*P
    dip_dur_time = dip_dur_phase*P

    for _t in t:
        ccf_resid_dip_model.append(
            fn_ccf_resid_RM_model_fixed_time(
                A_dip, vsini, t0_dip, dip_dur_time, _t, P, dv
            )
        )
    for _t in t_obs:
        ccf_resid_dip_obs.append(
            fn_ccf_resid_RM_model_fixed_time(
                A_dip, vsini, t0_dip, dip_dur_time, _t, P, dv
            )
        )

    ccf_resid_dip_model = np.array(ccf_resid_dip_model)
    ccf_resid_dip_obs = np.array(ccf_resid_dip_obs)

    # add in the residual to the full "dips" model
    ccf_resid_dips_model += ccf_resid_dip_model
    ccf_resid_dips_obs += ccf_resid_dip_obs


ccf_resid_model = ccf_resid_spot_model + ccf_resid_dips_model
np.random.seed(42)
gauss_err = np.random.normal(
    loc=0, scale=5e-3, size=ccf_resid_spot_obs.shape
)

red_err = (
    2e-3*(-t_obs[:,None]-t0)
    +
    1e-3*np.sin(3+1.1*24*t_obs[:,None])*(np.abs(velocity_grid[None,:]))**0.2
)

ccf_resid_obs = (
    ccf_resid_spot_obs + ccf_resid_dips_obs + gauss_err + red_err
)

#
# plot the ccf residuals
#

t_fn = lambda x: x*24

for plotid, _t, ccf_resid in zip(
    ['model', 'synthdata'],
    [t, t_obs],
    [ccf_resid_model, ccf_resid_obs]
):

    plt.close("all")
    set_style("science")
    factor=0.8
    fig, ax = plt.subplots(figsize=(factor*3,factor*4), layout='constrained')

    vscale = np.max([np.abs(ccf_resid.min()), np.abs(ccf_resid.max())])
    vmin = -vscale
    vmax = vscale

    c = ax.pcolor(velocity_grid,
                  t_fn(_t),
                  ccf_resid,
                  cmap='Spectral', vmin=vmin, vmax=vmax,
                  shading='auto')

    #ax.vlines(
    #    [-vsini, +vsini], t_fn(_t.min()), t_fn(_t.max()), zorder=10, lw=1, colors='k'
    #)

    colors = ['C0', 'C3', 'C2']
    for ix, (dip_ix, dip_vals) in enumerate(dip_dict.items()):
        _c = colors[ix]

        ax.vlines(
            2.5*vsini,
            t_fn(dip_dict[dip_ix]['t_ing']),
            t_fn(dip_dict[dip_ix]['t_egr']),
            #[velocity_grid.min(), velocity_grid.max()],
            zorder=10,
            lw=1,
            colors=_c
        )

        #ax.hlines(
        #    [t_fn(dip_dict[dip_ix]['t_ing']), t_fn(dip_dict[dip_ix]['t_egr'])],
        #    velocity_grid.min(),
        #    velocity_grid.max(),
        #    zorder=10,
        #    lw=1,
        #    colors=_c
        #)

    cb0 = fig.colorbar(c, ax=ax, extend='both',
                       location='top', shrink=0.5, pad=0.01,
                       ticks=[-0.03, 0, 0.03])
    cb0.set_label(f"{star_id} Synthetic Line Profile Residuals", labelpad=5)
    cb0.ax.set_xticklabels(['-3%', '0', '+3%'])

    ax.update({
        'xlabel':'Velocity [km/s]', 'ylabel':'Time [hours]',
        'ylim':[0, t_fn(_t.max())]
    })

    outpath = join(
        outdir, f"20230331_{star_id}_DT_simulation_ccf_resids_{plotid}.png"
    )
    savefig(fig, outpath, writepdf=0)
