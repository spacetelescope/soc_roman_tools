# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os
import copy
from warnings import warn
import numpy as np
import synphot as syn
from synphot.models import Empirical1D, ConstFlux1D
import astropy.units as u

from . import config as cf
from .constants import PANDEIA_WAVEUNITS, PANDEIA_FLUXUNITS, PANDEIA_ATOL, PANDEIA_RTOL
from .config import DefaultConfig
from .io_utils import ref_data_column
from .utils import spectrum_resample
from .custom_exceptions import EngineInputError, DataError, SynphotError, WavesetMismatch
import six

default_refdata_directory = cf.default_refdata_directory


class Background(cf.DefaultConfig):

    """
    This class provides the background flux.  The background can be computed from the engine's
    notional reference data or passed as a spectrum of a form [wave, sb] where wave is an array of
    wavelengths in microns and sb is an array of surface brightnesses in MJy/sr.

    Parameters
    ----------
    observation: Observation instance
        Contains needed instrument and background configuration data.
    config: dict
        Extra configuration information
    webapp: bool
        Enable strict API checking (not yet implemented here)

    Attributes
    ----------
    ref_dir: string
        Directory containing data for notional background model
    bg_level: string or list-like or array-like
        Input background level. Can be string for configuring notional model or lists/arrays for user-input background.
    instrument: Instrument instance
        Instrument configuration data is required to determine pixel scale and telescope/instrument contributions
        to the background model.
    aperture_pars: dict
        Instrument configuration data containing pixel scale required to calculate per-pixel background flux
    bg_spec: synphot.SourceSpectrum
        Original full-length, full-resolution background spectrum with wavelengths in microns and surface brightness in MJy/sr.
    bg: synphot.SourceSpectrum
        Current, possibly resampled background spectrum with wavelengths in microns and surface brightness in MJy/sr.

    Properties
    ----------
    wave: 1D np.ndarray
        Binned/trimmed wavelength set
    MJy_sr: 1D np.ndarray
        Binned/trimmed background spectrum in units of MJy/sr
    mjy_pix: 1D np.ndarray
        Convert self.MJy_sr to units of mJy (milli-janskys) per pixel
    flambda_pix: 1D np.ndarray
        Convert self.MJy_sr to units of F_lambda per pixel, photons/cm^2/s/micron/pixel.
    """

    def __init__(self, observation, config={}, webapp=False, **kwargs):
        self.webapp = webapp
        self.telescope = observation.instrument.telescope.tel_name

        DefaultConfig.__init__(self, config=config, webapp=webapp, **kwargs)

        self.ref_dir = os.path.join(default_refdata_directory, self.telescope, 'background')
        self.bg_location = observation.background
        self.instrument = observation.instrument
        self.aperture_pars = self.instrument.get_aperture_pars()
        # if the background spectrum is provided, the 0th index is wavelength and the 1st is the surface brightness
        if isinstance(self.bg_location, (list, tuple, np.ndarray)):
            try:
                # wavelength set for background that will get resampled to match wavelength set for a scene
                wave = np.array(self.bg_location[0], dtype=np.float64) << PANDEIA_WAVEUNITS
                # background spectrum in standard units of MJy/sr. use properties mjy_pix and flambda_pix to convert units.
                MJy_sr = np.array(self.bg_location[1], dtype=np.float64) << u.MJy
                # we store the original bg spectrum in a synphot spectrum so that resampling operations
                # are done on the original data.
                self.bg_spec = syn.spectrum.SourceSpectrum(Empirical1D, points=wave, lookup_table=MJy_sr)
            except Exception as e:
                msg = "Malformed input background spectrum: %s (%s)" % (repr(self.bg_location), e)
                raise EngineInputError(value=msg)
        elif isinstance(self.bg_location, (str, six.text_type)):
            if self.bg_location == "none":
                pass
            else:
                self.bg_level = observation.background_level
            self.get_canned_background()
        else:
            msg = "Must provide a valid background spectrum or background signifier string: %s" % repr(self.bg_location)
            raise EngineInputError(value=msg)
        self.bg = copy.deepcopy(self.bg_spec)

    def _get_config(self):
        """
        Read default configuration from JSON
        defaults are grouped based on telescope, get the background for the specified telescope

        Returns
        -------
        config: dict
            All desired class attributes should have defaults defined in the config file
        """
        # make use of the base class method to get the config
        config = super()._get_config()
        return config[self.telescope]

    def resample(self, wavelengths):
        """
        Re-bin background onto a new set of wavelengths.

        This method rewrites bg, always starting from the original unmodified bg_spec

        Parameters
        ----------
        wavelengths: 1D numpy array
            Array of new wavelengths to sample spectrum onto
        """
        # Yes, spectrum_resample expects input in mJy and outputs in mJy, but given that there's a constant
        # 1e9 factor also in fnu space between them, the results should be valid as MJy.
        MJy_sr = spectrum_resample(self.bg_spec(self.bg_spec.waveset, u.MJy).to_value(u.MJy), \
                                   self.bg_spec.waveset.to_value(PANDEIA_WAVEUNITS), wavelengths)
        self.bg = syn.spectrum.SourceSpectrum(Empirical1D, points=wavelengths << PANDEIA_WAVEUNITS, lookup_table=MJy_sr << u.MJy)

    def get_canned_background(self):
        """
        Read reference data for pandeia's canned background models and configure self.wave, self.bg_location,
        self.bg_level, and self.bg_spec to use it.
        """

        if self.bg_location == 'none' or self.bg_level == 'none':
            bg_wave = np.array([0.01, 30.0]) << PANDEIA_WAVEUNITS
            bg_sb = np.array([0.0, 0.0]) << PANDEIA_FLUXUNITS
            try:
                # Using ConstFlux1D would work, but then the SourceSpectrum has no waveset, and the rest of the functions expect one.
                tot_bg = syn.spectrum.SourceSpectrum(Empirical1D, points=bg_wave, lookup_table=bg_sb)
            except Exception as e:
                msg = "Problem creating background spectrum with Synphot."
                if self.webapp:
                    msg += "(%s)" % type(e)
                else:
                    msg += repr(e)
                raise SynphotError(value=msg)
        else:
            bg_file = os.path.join(default_refdata_directory, self.telescope, "background", "{0}_{1}.fits".format(self.bg_location,self.bg_level))
            try:
                bg_wave = ref_data_column(bg_file, colname='wavelength') << PANDEIA_WAVEUNITS
                bg = ref_data_column(bg_file, colname='background') << u.MJy #/u.sr, but that's not valid in spectra
                tot_bg = syn.spectrum.SourceSpectrum(Empirical1D, points=bg_wave, lookup_table=bg)
            except Exception as e:
                msg = "Problem loading background model. "
                if self.webapp:
                    msg += "(%s)" % type(e)
                else:
                    msg += repr(e)
                raise DataError(value=msg)

        self.bg_spec = tot_bg

    @property
    def wave(self):
        return self.bg.waveset.to_value(PANDEIA_WAVEUNITS)

    @property
    def MJy_sr(self):
        return self.bg(self.bg.waveset, flux_unit=u.MJy).to_value(u.MJy)


    @property
    def mjy_pix(self):
        """
        Convert background surface brightness in MJy/sr to F_nu per pixel in
        mJy/pixel (self.wave assumed to be in microns)

        Returns
        -------
        i_pix: 1D np.ndarray
            Flux per pixel in mJy/pixel
        """
        i_nu = self.bg(self.bg.waveset, flux_unit=PANDEIA_FLUXUNITS).to_value(PANDEIA_FLUXUNITS) * (u.arcsec * u.arcsec).to(u.sr)

        i_nu_pix = i_nu * self.aperture_pars['pix'] ** 2
        return i_nu_pix

    def trim(self, wmin, wmax):
        """
        Trim spectrum to wavelength range specified by wmin and wmax.
        This function does not use Astropy Quantities (except at the beginning and end)
        Like resample, it resamples off of the original background wavelength array

        Parameters
        ----------
        wmin: float
            Minimum wavelength (microns)
        wmax: float
            Maximum wavelength (microns)
        """
        wave = self.bg_spec.waveset.to_value(PANDEIA_WAVEUNITS)
        # Convert to MJy (synphot), output in units of MJy (astropy).
        flux = self.bg_spec(self.bg_spec.waveset, flux_unit=u.MJy).to_value(u.MJy)
        valid_subs = np.where((wave >= wmin) & (wave <= wmax))
        if valid_subs is None or len(valid_subs[0]) == 0:
            # corner case: We have a spectrum that covers the wavelength range 
            # of the instrument setup, but it's so sparsely specified that no 
            # points fall within the range.
            if (wave[0] <= wmin) and (wave[-1] >= wmax):
                wave = np.asarray([wmin,wmax])
                valid_subs = slice(2)
                flux = self.bg_spec(wave << PANDEIA_WAVEUNITS, flux_unit=u.MJy).to_value(u.MJy)
            else:
                msg = "Spectrum for source does not have any overlap with instrument wavelength range."
                raise WavesetMismatch(value=msg)

        trim_wave = []
        # make sure trimmed wavelengths include wmin and wmax
        if not np.isclose(wave[valid_subs][0], wmin, atol=PANDEIA_ATOL, rtol=PANDEIA_RTOL):
            trim_wave.append(wmin)

        trim_wave.extend(wave[valid_subs])

        if not np.isclose(wave[valid_subs][-1], wmax, atol=PANDEIA_ATOL, rtol=PANDEIA_RTOL):
            trim_wave.append(wmax)

        trim_wave = np.array(trim_wave)

        # this will use synphot to properly sample the endpoints
        # It requires unitless data
        trim_flux = spectrum_resample(flux, wave, trim_wave)

        self.nw = wave.size
        # re-add units
        self.bg = syn.spectrum.SourceSpectrum(Empirical1D, points=trim_wave << PANDEIA_WAVEUNITS, lookup_table=trim_flux << u.MJy)

    @property
    def flambda_pix(self):
        """
        Convert background surface brightness in MJy/sr to F_lambda per pixel in
        photons/cm^2/s/micron/pixel (self.wave assumed to be in microns)

        Returns
        -------
        i_pix: 1D np.ndarray
            Flux per pixel in photons/s/cm^2/pixel
        """
        i_lambda = 1.50916e9 / self.wave * self.bg(self.bg.waveset, flux_unit=syn.units.FLAM).value << (u.arcsec * u.arcsec).to(u.sr)
        i_lambda_pix = i_lambda * self.aperture_pars['pix'] ** 2
        return i_lambda_pix
