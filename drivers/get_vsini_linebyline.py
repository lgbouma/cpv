import os
from numpy import array as nparr
from scipy.interpolate import interp1d
from astropy.io import fits
from glob import glob
from os.path import join
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from cdips_followup.spectools import get_naive_rv
from complexrotators.paths import DATADIR, RESULTSDIR
from astropy import units as u
from aesthetic.plot import set_style, savefig

outdir = join(RESULTSDIR, 'HIRES_vsini_linebyline')
if not os.path.exists(outdir): os.mkdir(outdir)

def get_vsini_linebyline():

    starid = 'TIC141146667'
    VBROAD = 130
    fitsdir = join(DATADIR, 'spectra/HIRES/TIC141146667_DEBLAZED')
    run_in_parallel = 1

    for chip, order in zip(
        ['i'],
        [8]
    ):

        fitspaths = np.sort(glob(join(fitsdir, f"{chip}j*fits")))

        RUNDICT = {}
        mjds = []
        for fitspath in fitspaths:

            hdul = fits.open(fitspath)
            hdr = hdul[0].header
            mjd = hdr['MJD']
            mjds.append(mjd)
            hdul.close()

            ind = os.path.basename(fitspath).rstrip(".fits").lstrip(f"{chip}")
            key = f'TIC141146667_{ind}'
            # file name, grid teff, %2f grid logg, expected RV (if there is one)
            val = [os.path.basename(fitspath), 3000, 4.50, 0]
            RUNDICT[key] = val


        for starname,v in RUNDICT.items():

            fitsname, teff, logg, rv_expected = v

            dirstarname = starname
            if "_" in starname:
                dirstarname = starname.split("_")[0]

            spectrum_path = join(fitsdir, fitsname)

            localdir = join(os.path.expanduser('~'), 'local')
            synth_path0 = join(
                localdir,
                f'synthetic_spectra/PHOENIX_MedRes/'
                f'lte0{teff:d}-{logg:.2f}-0.0'
                f'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
            )

            # acquire RV and low-res template
            df, wav, flx, swav_lores, sflx_lores = get_naive_rv(
                spectrum_path, synth_path0, outdir, chip, make_plot=1,
                run_in_parallel=run_in_parallel, vbroad=VBROAD,
                specific_orders=[order], return_spectra=1,
                overwrite=1
            )

            synth_path1 = join(
                localdir,
                f'synthetic_spectra/PHOENIX_HiRes/'
                f'lte0{teff:d}-{logg:.2f}-0.0'
                f'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
            )

            from cdips_followup.spectools import (
                get_synth_spectrum, degrade_spectrum, broaden_spectrum
            )
            syn_flx, syn_wav = get_synth_spectrum(synth_path1)

            min_wav, max_wav = np.nanmin(wav), np.nanmax(wav)
            wav_0 = np.nanmean(wav)*u.AA
            syn_sel = (syn_wav > min_wav) & (syn_wav < max_wav)
            sflx = syn_flx[syn_sel] / np.nanmedian(syn_flx[syn_sel])
            swav = syn_wav[syn_sel]
            fudge = 0.9
            dlam = fudge * (
                float(df['rv_chisq_kms']) * wav_0 /
                (299792.458)
            ).value

            # the template / synthetic spectrum observed at HIRES
            # instrumental resolution
            sflx_hiresR = degrade_spectrum(swav, sflx, R=45000)

            # the template / synthetic spectrum observed at a range of
            # vsini's
            vsinis = np.arange(130, 170+10, 10)
            sflx_vsinis = []
            for vsini in vsinis:
                epsilon = 1.5 # linear limb darkening coeff
                # broaden
                bsflx_vsini = broaden_spectrum(
                    swav, sflx, vsini=vsini, epsilon=epsilon
                )
                # then interpolate down to observed wavelength grid
                interpolator = interp1d(swav, bsflx_vsini, kind="quadratic",
                                        fill_value="extrapolate")
                sflx_vsini = interpolator(wav)
                sflx_vsinis.append(sflx_vsini)

            N = len(vsinis)
            colors = [plt.cm.viridis_r(i) for i in np.linspace(0, 0.8, N)]

            ##########
            dy = 0.

            set_style('science')
            plt.close('all')

            f = 1
            fig, ax = plt.subplots(figsize=(f*4,f*3))

            # broaden the data spectrum by much less than the target
            # wavelength for visual comparison
            bflx = broaden_spectrum(wav, flx, vsini=10,
                                    epsilon=epsilon)

            norm_wv0 = 7701 # continuum normalization is wild
            norm_wv1 = 7712
            ind_lo = np.argmin(np.abs((wav - dlam) - norm_wv0))
            ind_hi = np.argmin(np.abs((wav - dlam) - norm_wv1))
            b_lo = bflx[ind_lo]
            b_hi = bflx[ind_hi]
            A = b_hi / b_lo

            ax.plot(wav-dlam, bflx, c='k', lw=1.5, zorder=99)

            for ix, (sflx_vsini, vsini, c) in enumerate(
                zip(sflx_vsinis, vsinis, colors)
            ):
                w0, w1 = 7700, 7713
                sel = (wav > w0) & (wav < w1)

                # Set yscale as line depth vs "continuum" just outside
                # and offset in y to match continuum
                # NOTE: I tried the usual idea of a least-squares
                # solution and it did not work.
                s_ind_lo = np.argmin(np.abs(wav - norm_wv0))
                s_ind_hi = np.argmin(np.abs(wav - norm_wv1))

                s_lo = sflx_vsini[s_ind_lo]
                s_hi = sflx_vsini[s_ind_hi]
                B = s_hi / s_lo
                factor = A/B

                offset = (factor*sflx_vsini)[s_ind_lo] - bflx[ind_lo]

                print(f'{ix}: factor={factor:.3f}, offset={offset:.3f}')

                ax.plot(wav, factor*sflx_vsini - offset, c=c, lw=1,
                        label=f'{vsini}', zorder=ix)
            #TODO FIXME OK GREAT; ONCE THIS IS DONE, DO CHISQ  OR LSQ
            #ERR FOR EACH OF THESE BTWN YOUR WAVELENGTH REGIONS...
            #TODO FIXME OK GREAT; ONCE THIS IS DONE, DO CHISQ  OR LSQ
            #ERR FOR EACH OF THESE BTWN YOUR WAVELENGTH REGIONS...
            #TODO FIXME OK GREAT; ONCE THIS IS DONE, DO CHISQ  OR LSQ
            #ERR FOR EACH OF THESE BTWN YOUR WAVELENGTH REGIONS...
            #TODO FIXME OK GREAT; ONCE THIS IS DONE, DO CHISQ  OR LSQ
            #ERR FOR EACH OF THESE BTWN YOUR WAVELENGTH REGIONS...

            ax.legend(fontsize='xx-small')
            ax.update({'xlabel':'λ', 'ylabel':'f'})

            ax.set_xlim([7695.5, 7713])
            ax.set_ylim([0.45, 1.25])

            tname = os.path.basename(spectrum_path).replace(".fits", "")
            sname = os.path.basename(synth_path1).replace(".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits", "")
            odir = os.path.join(outdir, tname+"_v_"+sname)
            if not os.path.exists(odir): os.mkdir(odir)

            savpath = join(
                odir, f'vsini_ord{str(order).zfill(2)}_{tname}.png'
            )
            savefig(fig, savpath)


if __name__ == "__main__":
    get_vsini_linebyline()


