
from pandeia.engine.helpers.background.jwst.sl_accessor_funcs import get_dateless_bg
from pandeia.engine.helpers.background.jwst.straylight_accessor import get_stray_light_bg
from pandeia.engine.helpers.background.jwst.thermal_accessor import get_thermal_bg
from pandeia.engine.helpers.background.jwst.bmg_accessor import get_in_field_bg
from pandeia.engine.helpers.background.hst.sky import Sky
from pandeia.engine.helpers.background.hst.sky_spectrum import SkySpectrum

from pandeia.engine.custom_exceptions import BMGError
from pandeia.engine.custom_exceptions import BackgroundError
from pandeia.engine.custom_exceptions import DatelessBGError
from pandeia.engine.custom_exceptions import StraylightPositionError
from pandeia.engine.custom_exceptions import StraylightDataError

import synphot as syn
import astropy.units as u
from synphot.models import ConstFlux1D, Empirical1D
from operator import add
import numpy as np
import os

'''
We have to resample the three background components because the infield waveset
is not the same as straylight and thermal.  So we normalize the waveset and then
resample each component.
'''
def bg_resample(merged, wave, flux):
    merged = np.array(merged)
    wave = np.array(wave)
    flux = np.array(flux)
    # we are NOT converting these spectra to microns and MJy (which should be MJy/sr anyway)
    # because that will trigger the internal unit conversion routines, and they have been found
    # to corrupt the spectra. Because we do not need to do any actual unit conversions, we can 
    # allow synphot to think these are its default units (Angstroms and PHOTLAM) and simply 
    # accept the rebinned values. (JETC-300) See similar code in engine/utils.py spectrum_resample
    spec = syn.spectrum.SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)
    filt = syn.spectrum.SpectralElement(ConstFlux1D, amplitude=1)
    obs = syn.observation.Observation(spec, filt, binset=merged, force='taper')
    flux = obs.binflux.value
    return list(flux)


def call_butler(ra, dec, date, level, ra_dec_str, date_str):

    if date is not None:
        # dated

        # stray light
        try:
            sl_wave, sl_bg = get_stray_light_bg(ra, dec, date, ra_dec_str, date_str)
        except (StraylightPositionError, StraylightDataError):
            raise
        except Exception as e:
            raise BackgroundError('Error calculating stray light background.', e)

        # thermal
        try:
            thermal_wave, thermal_bg = get_thermal_bg()
        except Exception as e:
            raise BackgroundError('Error calculating thermal background.', e)

        # in field
        try:
            if_wave, if_bg = get_in_field_bg(ra, dec, date, ra_dec_str, date_str)
        except BMGError as e:
            raise BackgroundError('BMGError calculating infield background.  This can happen if the BMG is unavailable, or because %s.' % e)
        except Exception as e:
            if "timed out" in str(e):
                raise BackgroundError("Error calculating infield background: dated background is currently unavailable.  Try again later or switch to low/medium/high")
            else:
                raise BackgroundError('Error calculating infield background: %s' % e)

        # merge wavelengths
        # wave = syn.units.merge_wavelengths(thermal_wave, if_wave)
        wave = if_wave

    else:

        # get all components via special get_dateless_bg route if no date given
        try:
            unused_iday, wave, if_bg, sl_bg, thermal_wave, thermal_bg = get_dateless_bg(ra, dec, level)
            # the get_dateless_bg() API returns iday (int day of yr) which we don't currently use
        except (DatelessBGError, StraylightPositionError, StraylightDataError):
            raise
        except Exception as e:
            raise BackgroundError('Error calculating dateless background.', e)
        # wave stands for both if_wave (infield) and sl_wave (straylight)

    # resample the flux for each onto the new set
    if date is not None:
        sl_bg = bg_resample(wave, sl_wave, sl_bg)
    thermal_bg = bg_resample(wave, thermal_wave, thermal_bg)
    # if_bg = bg_resample(wave, if_wave, if_bg)

    # then sum the backgrounds
    combined_bg = list(map(add, thermal_bg, if_bg))
    combined_bg = list(map(add, combined_bg, sl_bg))

    # this will be written to a .npz file
    data_to_save = dict(
        straylight=sl_bg,
        thermal=thermal_bg,
        infield=if_bg,
        background=combined_bg,
        wavelength=wave
    )

    return [list(wave), combined_bg], data_to_save


def get_jwst_background(background):
    data_to_save = None

    # dated
    if background['bg_type'] == 'dated':
        ra = background['ra']
        dec = background['dec']
        ra_dec_str = background['ra_dec_str']
        date = background['date']
        date_str = background['date_str']
        background, data_to_save = call_butler(ra, dec, date, None, ra_dec_str, date_str)

    # dateless
    elif background['bg_type'] is not None and background['bg_type'].lower() != 'none':
        ra = background['ra']
        dec = background['dec']
        ra_dec_str = background['ra_dec_str']
        level = background['bg_type'].upper()[0]  # dateless code expects level of: L,M,H
        background, data_to_save = call_butler(ra, dec, None, level, ra_dec_str, "dateless")

    else:
        # Handle "none", but note that this is how we'd send positionless/dateless
        background = background['bg_type']

    return background, data_to_save


def get_hst_background(bg_input):
    data_to_save = None
    background = np.array([0.01, 30.0]), np.array([0.0, 0.0])

    bg = Sky(bg_input)
    bg.make_synexpr()

    if bg.expr is not None:
        bg_spectrum = SkySpectrum(bg.expr)
        # convert from angstroms and flam/arcsec^2 to microns and MJy/sr
        wavelength = bg_spectrum.wave_in_microns
        flux = bg_spectrum.flux_in_MJy_per_sr
        background = wavelength, flux

        # this will be written to a .fits file
        data_to_save = dict(
            background=flux,
            wavelength=wavelength
        )

    return background, data_to_save


def get_background(background):
    data_to_save = None
    telescope = os.environ.get("ETC_TELESCOPE", "jwst").lower()

    if telescope == 'hst':
        background, data_to_save = get_hst_background(background)
    else:
        background, data_to_save = get_jwst_background(background)

    return background, data_to_save

