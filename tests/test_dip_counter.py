import os
from os.path import join
from complexrotators.paths import LOCALDIR
from complexrotators.plotting import plot_phased_light_curve
from complexrotators.getters import _get_lcpaths_given_ticid
from complexrotators.lcprocessing import (
    cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
)

def test_dip_counter(ticid):

    lcpaths = _get_lcpaths_given_ticid(ticid)
    cachedir = join(LOCALDIR, "cpv_finding")
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    for lcpath in lcpaths:

        # get the relevant light curve data
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        # get period, t0, and periodogram (PDM or LombScargle)
        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        r = count_phased_local_minima(
            d['times'], d['fluxs'], d['t0'], d['period'],
            binsize_phase_units=0.005, prominence=1e-3, width=3
        )

        import IPython; IPython.embed()


def test_dip_counter_all_stars():

    ticids = [
    #"201789285",
    #"311092148",
    #"332517282",
    #"405910546",
    #"142173958",
    #"300651846",
    #"408188366",
    #"146539195",
    #"177309964",
    #"425933644",
    #"206544316",
    #"224283342",
    "245902096",
    #"150068381",
    #"177309964",
    #"118769116",
    #"245868207",
    #"245874053",
    #"59129133"
    ]

    for ticid in ticids:
        test_dip_counter(ticid)

if __name__ == "__main__":
    test_dip_counter_all_stars()
