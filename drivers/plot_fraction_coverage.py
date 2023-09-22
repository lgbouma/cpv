"""
Make plot of completeness fraction for the 2-minute data between Sectors 1-55.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt

def make_figure(factor=2):

    # S1-S58 2-min data
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

    N_TESS_2min = len(np.unique(stdf.dr2_source_id))

    # Gaia 100pc
    do_100pc = 0
    do_150pc = 1
    if do_100pc:
        from rudolf.helpers import get_gaia_catalog_of_nearby_stars
        gdf = get_gaia_catalog_of_nearby_stars()
    elif do_150pc:
        # Gaia DR2 <150pc
        # select
        # g.source_id, g.ra, g.dec, g.parallax, g.parallax_over_error, g.pmra,
        # g.pmdec, g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
        # g.phot_g_mean_flux_over_error, g.phot_bp_mean_flux_over_error,
        # g.phot_rp_mean_flux_over_error, g.bp_rp, g.g_rp, g.radial_velocity, g.l,
        # g.b
        # from gaiadr2.gaia_source import as g
        # where g.parallax > 6.6666666666
        # and g.bp_rp > 1.5
        gdf = pd.read_csv('/Users/luke/local/complexrotators/lt150pc_mkdwarf-result.csv')


    Tmag_pred = (
            gdf['phot_g_mean_mag']
            - 0.00522555 * (gdf['phot_bp_mean_mag'] - gdf['phot_rp_mean_mag'])**3
            + 0.0891337 * (gdf['phot_bp_mean_mag'] - gdf['phot_rp_mean_mag'])**2
            - 0.633923 * (gdf['phot_bp_mean_mag'] - gdf['phot_rp_mean_mag'])
            + 0.0324473
    )
    gdf['TESSMAG'] = Tmag_pred

    gdf['M_G'] = gdf['phot_g_mean_mag'] + 5*np.log10(gdf['parallax']/1e3) + 5

    sel = (
        (gdf.phot_bp_mean_mag - gdf.phot_rp_mean_mag > 1.5)
        &
        (gdf.TESSMAG < 16)
        &
        (gdf.parallax > 1e3*(1/150))
        &
        (gdf.M_G > 4)
    )
    sgdf = gdf[sel]

    N_Gaia = len(np.unique(sgdf.source_id))

    frac = 100 * N_TESS_2min / N_Gaia

    print(f'{frac:.1f}% ({N_TESS_2min} / {N_Gaia}) of K- and M-dwarfs with d<150pc and T<16 received 2min '
          'observations acrosrs Sectors 1-55')

    delta_dist = 10
    dist_bins = np.arange(0, 160, delta_dist)
    dist_mids = dist_bins[0:-1] + delta_dist/2

    fracs = []
    nums, denoms = [], []

    for ix, dist_bin in enumerate(dist_bins):

        if dist_bin == max(dist_bins):
            continue

        dist_low = dist_bin
        dist_high = dist_bins[ix+1]
        print(ix, dist_low, dist_high)

        if ix != 0 :
            get_sel = (
                lambda x:
                (x.phot_bp_mean_mag - x.phot_rp_mean_mag > 1.5)
                &
                (x.TESSMAG < 16)
                &
                (x.M_G > 4)
                &
                (x.parallax > 1e3*(1/dist_high))
                &
                (x.parallax <= 1e3*(1/dist_low))
            )
        else:
            get_sel = (
                lambda x:
                (x.phot_bp_mean_mag - x.phot_rp_mean_mag > 1.5)
                &
                (x.TESSMAG < 16)
                &
                (x.M_G > 4)
                &
                (x.parallax > 1e3*(1/dist_high))
            )

        stdf = tdf[get_sel(tdf)]
        sgdf = gdf[get_sel(gdf)]

        N_TESS_2min = len(np.unique(stdf.dr2_source_id))
        N_Gaia = len(np.unique(sgdf.source_id))

        nums.append(N_TESS_2min)
        denoms.append(N_Gaia)
        fracs.append(N_TESS_2min/N_Gaia)

    outdf = pd.DataFrame({
        'dist_mids_pc': dist_mids,
        'N_TESS_2min': nums,
        'N_Gaia': denoms,
        'frac': fracs
    })
    if factor == 2:
        outdf.to_csv('../results/fraction_coverage/fraction_coverage.csv', index=False)

    # make fraction vs distance plot
    from aesthetic.plot import set_style, savefig

    plt.close("all")
    set_style("science")
    fig, ax = plt.subplots(figsize=(factor*2, factor*2))
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

    ax2 = ax.twinx()
    ax2.plot(outdf.dist_mids_pc, outdf.N_Gaia, c='C1', marker='o', lw=2, zorder=-3,
             ms=factor*4.5, ls=':')
    ax2.set_yscale('log')

    ax.plot(outdf.dist_mids_pc, outdf.frac, c='C0', marker='P', lw=2, zorder=4,
            ms=factor*4.5)

    ax.update({
        'xlabel': 'Distance [pc]',
        'ylim': [0, 1]
    })
    ax.tick_params(axis='y', labelcolor='C0')
    ax.set_ylabel('Fraction with 2-min', color='C0')

    ax2.set_ylabel(
        'Stars per 10pc shell',
        rotation=270,
        color='C1',
        labelpad=10
    )
    ax2.tick_params(axis='y', labelcolor='C1')

    ax.set_title(
        '$T$<16, $G_{\mathrm{BP}}$-$G_{\mathrm{RP}}$>$1.5$, $M_{\mathrm{G}}$>$4$',
        pad=10
    )

    # default for manuscript
    if factor == 2:
        savefig(fig, "../results/fraction_coverage/fraction_coverage.png")
    else:
        savefig(
            fig,
            f"../results/fraction_coverage/fraction_coverage_factor{factor}.png"
        )

if __name__ == "__main__":
    make_figure(factor=1.5)
    make_figure(factor=2)
