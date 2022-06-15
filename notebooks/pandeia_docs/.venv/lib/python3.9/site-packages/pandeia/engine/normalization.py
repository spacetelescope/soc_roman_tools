# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os
import numpy as np
import astropy.units as u

# suppress stsynphot's unnecessary warning (we'll supply our own Vega spectrum)
import warnings
from astropy.utils.exceptions import AstropyUserWarning
# suppress stsynphot's unnecessary warning
warnings.filterwarnings("ignore", "Failed to load Vega spectrum from", AstropyUserWarning)

import synphot as syn
import stsynphot as stsyn
from synphot.models import Empirical1D

from .config import DefaultConfig
from .custom_exceptions import DataConfigurationError, EngineInputError, SynphotError
from .pandeia_warnings import normalization_warning_messages as warning_messages
from .utils import get_key_list, get_dict_from_keys, recursive_subclasses, merge_data
from .instrument_factory import InstrumentFactory
from .constants import PANDEIA_WAVEUNITS, PANDEIA_FLUXUNITS, UNIT_MAP

from . import io_utils as io
from . import config as cf

default_refdata_directory = cf.default_refdata_directory

class Normalization(DefaultConfig):

    """
    Class for handling configuration and application of spectrum normalization

    Attributes
    ----------
    type: string
        Type of normalization to perform
    norm_fluxunit: string
        Unit for fluxes used in normalization
    norm_waveunit: string
        Unit for wavelength used for reference wavelength, 'norm_wave'
    norm_wave: float
        Reference wavelength in 'norm_waveunit' at which the spectrum will be scaled for type="at_lambda"
    norm_flux: float
        Reference flux at 'norm_wave' to which spectrum is scaled to match for type="at_lambda"
    bandpass: string
        Specifies how bandpass information is looked up or calculated for normalization
        types 'hst', 'jwst', and 'photsys'.
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        How a normalization is configured depends on the defaults, the per-type defaults, and any input parameters.

        Parameters
        ----------
        webapp: bool
            Determines whether strict API checks should be performed
        config: dict
            Extra configuration information in engine input API dict format
        **kwargs: keyword/value pairs
            Extra configuration information in kwargs format
        """
        self.webapp = webapp

        # send webapp=False here since we do the API checking further below
        DefaultConfig.__init__(self, webapp=False, config=config, **kwargs)

        all_config = merge_data(config, dict(**kwargs))

        if not hasattr(self, "type"):
            self.type = self.__class__.__name__.lower().replace('normalize', '')

        if self.type in self.types:
            type_conf = self.types[self.type]
            if "defaults" in type_conf:
                type_defaults = type_conf.pop("defaults")
            else:
                type_defaults = {}
            type_conf = merge_data(type_conf, type_defaults, config, dict(**kwargs))
            self.__dict__ = merge_data(self.__dict__, type_conf)

        if hasattr(self,'norm_waveunit'):
            try:
                # If it's not already a quantity, it will be turned into one
                self.norm_wave = syn.units.validate_quantity(self.norm_wave, UNIT_MAP[self.norm_waveunit])
            except Exception as err:
                msg = "Failed to convert reference wavelength units, {}, to pandeia wavelength units, {}: ".format(self.norm_waveunit, PANDEIA_WAVEUNITS, err)
                if self.webapp:
                    msg += "{}".format(type(err))
                else:
                    msg += repr(err)
                raise EngineInputError(value=msg)

        if hasattr(self,'norm_fluxunit'):
            try:
                # If it's not already a quantity, it will be turned into one
                self.norm_flux = syn.units.validate_quantity(self.norm_flux, UNIT_MAP[self.norm_fluxunit])
            except Exception as err:
                msg = "Failed to convert reference flux units, {}, to pandeia flux units, {}: ".format(self.norm_fluxunit, PANDEIA_FLUXUNITS)
                if self.webapp:
                    msg += "{}".format(type(err))
                else:
                    msg += repr(err)
                raise EngineInputError(value=msg)

        # do some sanity checks to make sure inputs make sense
        try:
            self._sanity_checks()
        except AttributeError as e:
            self.warnings['no_sanity_check'] = warning_messages['no_sanity_check'] % (self.__class__.__name__, e)

        # we need to run the API checks here after merging in the normalization type-specific defaults
        if webapp:
            self._api_checks(all_config)

    def _get_config(self):
        """
        Read default configuration from JSON

        Returns
        -------
        config: dict
            All desired class attributes should have defaults defined in the config file
        """
        # use this trick to key the configuration file name off the name of the instantiated subclass
        ref_dir = os.path.join(default_refdata_directory, "normalization")
        config = io.read_json(os.path.join(ref_dir, "config.json"), raise_except=True)

        # pop the defaults entry out
        if "defaults" in config:
            defaults = config.pop("defaults")
            config.update(defaults)
        else:
            msg = "No normalization defaults defined."
            raise DataConfigurationError(value=msg)

        return config

    def _sanity_checks(self):
        """
        Make sure norm_fluxunit is supported and type is a configured normalization type
        """
        if self.norm_fluxunit not in self.fluxunits:
            msg = "Unsupported flux unit %s for normalization type %s" % (self.norm_fluxunit, self.type)
            raise EngineInputError(value=msg)

        # Make sure the unit matches what norm_fluxunit is supposed to be
        if self.norm_flux.unit != UNIT_MAP[self.norm_fluxunit]:
            msg = "Unit mismatch: flux is {}, supplied unit is {}".format(self.norm_flux.unit, self.norm_fluxunit)
            raise EngineInputError(value=msg)

        if self.type not in self.types:
            msg = "No configuration data found for normalization type %s" % self.type
            raise DataConfigurationError(value=msg)

class NormalizeAtLambda(Normalization):

    """
    Subclass for handling normalization via a reference flux, self.norm_flux, at a reference wavelength, self.norm_wave.
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        Parameters
        ----------
        webapp: bool
            Determines whether strict API checks should be performed
        config: dict
            Extra configuration information in engine input API dict format
        **kwargs: keyword/value pairs
            Extra configuration information in kwargs format
        """
        Normalization.__init__(self, webapp=webapp, config=config, **kwargs)

        # check sanity of input units
        if self.norm_waveunit not in self.waveunits:
            msg = "Unsupported reference wavelength unit, %s. Supported units are: %s" % (self.norm_waveunit, repr(self.waveunits))
            raise EngineInputError(value=msg)

        if self.norm_wave.unit != UNIT_MAP[self.norm_waveunit]:
            msg = "Unit mismatch: wavelength is {}, supplied unit is {}".format(self.norm_wave.unit, self.norm_waveunit)
            raise EngineInputError(value=msg)

        if self.norm_fluxunit not in self.fluxunits:
            msg = "Unsupported flux unit, %s. Supported units are: %s" % (self.norm_fluxunit, repr(self.fluxunits))
            raise EngineInputError(value=msg)

        # use synphot to convert self.norm_flux to PANDEIA_FLUXUNITS, if necessary.
        try:
            self.norm_flux = syn.units.convert_flux(self.norm_wave, self.norm_flux, PANDEIA_FLUXUNITS)
        except Exception as e:
            msg = "Error using Synphot to convert norm_fluxunit to PANDEIA_FLUXUNITS. "
            if self.webapp:
                msg += "(%s)" % type(e)
            else:
                msg += repr(e)
            raise SynphotError(value=msg)

    def normalize(self, wave, flux):
        """
        Normalize spectrum by scaling it at self.norm_wave.

        Parameters
        ----------
        wave: 1D np.ndarray
            Wavelength vector in 'PANDEIA_WAVEUNITS'
        flux: 1D np.ndarray
            Flux vector in 'PANDEIA_FLUXUNITS' containing the spectrum

        Returns
        -------
        scaled_wave, scaled_flux: 1D np.ndarray
            Wavelength and scaled flux vectors in PANDEIA_WAVEUNITS and PANDEIA_FLUXUNITS
        """
        if self.norm_wave > wave.max() or self.norm_wave < wave.min():
            msg = "Specified normalization wavelength, {:.2f} microns, not within wavelength bounds of spectrum: ({:.2f}, {:.2f} " \
                "microns)".format(self.norm_wave.value, wave.min().to_value(PANDEIA_WAVEUNITS), wave.max().to_value(PANDEIA_WAVEUNITS))
            raise EngineInputError(value=msg)

        spec_flux = np.interp(self.norm_wave, wave, flux)
        if spec_flux == 0.0:
            key = 'normalized_to_zero_flux'
            self.warnings[key] = warning_messages[key]
            scale_factor = 1.0
        else:
            # This operation is NOT astropy unit aware, so we need to be sure both are in PANDEIA_FLUXUNITS.
            # They already should be, but this forces it.
            self.norm_flux = syn.units.validate_quantity(self.norm_flux, PANDEIA_FLUXUNITS)
            spec_flux = syn.units.validate_quantity(spec_flux, PANDEIA_FLUXUNITS)
            scale_factor = self.norm_flux / spec_flux
            normed_flux = np.interp(self.norm_wave, wave, flux/spec_flux)

        scaled_wave = np.copy(wave)
        scaled_flux = flux * scale_factor

        return scaled_wave, syn.units.convert_flux(scaled_wave, scaled_flux, PANDEIA_FLUXUNITS)


class NormalizePhotsys(Normalization):

    """
    Subclass for normalizing a spectrum to a reference flux/magnitude in a filter in a supported photometric system
    (e.g. Cousins, Bessell, SDSS)
    """

    def _get_bandpass(self, *args):
        """
        Parse a self.bandpass of the form <photsys>,<filter> and use the keys to look up a bandpass file
        in self.photsystems.  Use synphot to load the bandpass and return it.

        Returns
        -------
        bp: synphot.spectrum.SpectralElement
            Bandpass throughput data as returned by synphot.spectrum.SpectralElement.from_file() and converted to pandeia
            wavelength units
        """
        try:
            keys = get_key_list(self.bandpass, separator=',')  # should be [photsys, filter]
            bp_filename = get_dict_from_keys(self.bandpasses, keys)['filename']
            bp_path = os.path.join(os.environ.get("PYSYN_CDBS"), "comp", "nonhst", bp_filename)
            bp = syn.spectrum.SpectralElement.from_file(bp_path)
        except Exception as e:
            msg = "Error loading normalization bandpass: %s " % self.bandpass
            if self.webapp:
                msg += "(%s)" % type(e)
            else:
                msg += repr(e)
            raise DataConfigurationError(value=msg)
        return bp

    def normalize(self, wave, flux):
        """
        Normalize a spectrum to a bandpass fetched from self._get_bandpass()

        Parameters
        ----------
        wave: 1D np.ndarray
            Wavelength vector in pandeia wavelength units
        flux: 1D np.ndarray
            Flux vector in 'self.norm_fluxunit' containing the spectrum

        Returns
        -------
        scaled_wave, scaled_flux: 1D np.ndarrays
            Wavelength and scaled flux vectors in pandeia units, microns and mJy
        """

        sp = syn.spectrum.SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)
        bp = self._get_bandpass(wave)
        try:
            sp_rn = sp.normalize(self.norm_flux, band=bp, vegaspec=syn.spectrum.SourceSpectrum.from_vega())
        except Exception as e:
            msg = "Error using Synphot to perform bandpass normalization in %s" % str(bp)
            if self.webapp:
                msg += "(%s)" % type(e)
            else:
                msg += repr(e)
            raise SynphotError(value=msg)

        return sp_rn.waveset.to(PANDEIA_WAVEUNITS), sp_rn(sp_rn.waveset, flux_unit=PANDEIA_FLUXUNITS)


class NormalizeObsmode(NormalizePhotsys):

    """
    Subclass for normalizing a spectrum based on a synphot-compatible obsmode string via synphot.spectrum.SpectralElement()
    """

    def _get_bandpass(self, *args):
        """
        Wrap stsynphot.band() to generate a bandpass based on a valid obsmode specification

        Returns
        -------
        bp: synphot.SpectralElement
            synphot bandpass
        """
        if self.webapp and self.bandpass not in self.bandpasses:
            key = "unsupported_normalization_bandpass"
            self.warnings[key] = warning_messages[key] % self.bandpass

        try:
            bp = stsyn.band(str(self.bandpass))
        except Exception as e:
            msg = "Error using Synphot to load bandpass via ObsMode string, %s. " % self.bandpass
            if self.webapp:
                msg += "(%s)" % type(e)
            else:
                msg += repr(e)
            raise SynphotError(value=msg)
        return bp


class NormalizeHst(NormalizeObsmode):

    """
    Largely an alias for NormalizeObsmode since HST normalization bands will be specified as obsmode strings.
    """

    pass


class NormalizeJwst(NormalizePhotsys):

    """
    Subclass for normalizing a spectrum based on a JWST configuration
    """

    def _get_bandpass(self, wave):
        """
        Get JWST instrument, mode, filter from self.bandpass, use that to create a JWSTInstrument instance, and then
        query that for throughput vs wavelength.  Use synphot.SpectralElement to convert the results to
        synphot-compatible form.

        Parameters
        ----------
        wave: 1D np.ndarray
            Wavelength vector onto which JWST bandpass is interpolated

        Returns
        -------
        bp: synphot.SpectralElement
            synphot bandpass converted to pandeia wavelength units.
        """
        keys = get_key_list(self.bandpass, separator=',')
        if len(keys) != 3:
            msg = "JWST bandpass specification must be of the form <instrument>,<mode>,<filter>"
            raise EngineInputError(value=msg)

        instrument, mode, filt = keys
        config = {}
        config['instrument'] = {}
        config['instrument']['instrument'] = instrument
        config['instrument']['mode'] = mode
        config['instrument']['filter'] = filt
        inst = InstrumentFactory(config=config)
        thruput = inst.get_total_eff(wave.to_value(PANDEIA_WAVEUNITS))
        try:
            bp = syn.spectrum.SpectralElement(Empirical1D, points=wave, lookup_table=thruput)
        except Exception as e:
            msg = "Error creating Synphot bandpass from JWST bandpass specification, %s. " % self.bandpass
            if self.webapp:
                msg += "(%s)" % type(e)
            else:
                msg += repr(e)
            raise SynphotError(value=msg)
        return bp


class NormalizeNone(Normalization):

    """
    Subclass to make no normalization look like the other methods.
    """

    def normalize(self, wave, flux):
        """
        Implement this method to maintain consistent API and simply return what's given.

        Parameters
        ----------
        wave: 1D np.ndarray
            Wavelength vector in pandeia wavelength units
        flux: 1D np.ndarray
            Flux vector in 'self.norm_fluxunit' containing the spectrum

        Returns
        -------
        wave, flux: 1D np.ndarrays
            Return wave and flux without modification
        """
        return wave, syn.units.convert_flux(wave, flux, PANDEIA_FLUXUNITS)

    def _sanity_checks(self):
        """
        'type' is the only parameter this normalization type accepts. Should never get here
        if NormalizationFactory is used...
        """
        if self.type not in self.types:
            msg = "No configuration data found for normalization type %s" % self.type
            raise DataConfigurationError(value=msg)


def NormalizationFactory(webapp=False, config={}, **kwargs):

    """
    Function to take normalization configuration data and build/return a configured instance of the appropriate
    subclass of Normalization.

    Parameters
    ----------
    webapp: bool
        Switch to toggle strict API checking
    config: dict
        Configuration data in engine API dict format
    **kwargs: keyword/value pairs
        Additional configuration data
    """
    all_config = merge_data(config, dict(**kwargs))
    if 'type' in all_config:
        method = all_config['type'].replace('_', '')
    else:
        msg = "Must specify type of normalization to perform."
        raise EngineInputError(value=msg)

    types = recursive_subclasses(Normalization)
    methods = [t.__name__.lower().replace('normalize', '') for t in types]
    type_map = dict(list(zip(methods, types)))

    if method not in methods:
        msg = "Unsupported or not yet implemented normalization method: %s" % method
        raise EngineInputError(value=msg)
    else:
        cls = type_map[method](webapp=webapp, config=config, **kwargs)
        return cls
