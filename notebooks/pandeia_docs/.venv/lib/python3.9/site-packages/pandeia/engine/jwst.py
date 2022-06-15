# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os

from astropy.io import fits

import numpy.ma as ma
import numpy as np

from .telescope import Telescope
from .instrument import Instrument
from .io_utils import read_json
from .custom_exceptions import EngineInputError, RangeError, DataError

class JWST(Telescope):

    """
    This is currently a dummy class that is used for configuration file discovery. Could eventually
    contain JWST-specific methods.
    """
    pass


class JWSTInstrument(Instrument):

    """
    Generic JWST Instrument class
    """

    def __init__(self, mode=None, config={}, webapp=False, **kwargs):
        telescope = JWST()
        # these are the required sections and need to be passed via API in webapp mode
        self.instrument_pars = {}
        self.instrument_pars['detector'] = ["nexp", "ngroup", "nint", "readout_pattern", "subarray"]
        self.instrument_pars['instrument'] = ["aperture", "disperser", "filter", "instrument", "mode"]
        self.api_parameters = list(self.instrument_pars.keys())

        # these are required for calculation, but ok to live with config file defaults
        if hasattr(self, 'api_ignore'):
            self.api_ignore.extend(['dynamic_scene', 'max_scene_size', 'scene_size'])
        else:
            self.api_ignore = ['dynamic_scene', 'max_scene_size', 'scene_size']

        Instrument.__init__(self, telescope=telescope, mode=mode, config=config, webapp=webapp, **kwargs)

    def read_detector(self):
        """Read in the detector keyword from the aperture parameter in the config json file.
           Put the detector keyword in self.instrument['detector']."""
        if "slit" in self.aperture_config[self.get_aperture()]:
            self.instrument['slit'] = self.aperture_config[self.get_aperture()]['slit']
        self.instrument['detector'] = self.aperture_config[self.get_aperture()]['detector']

class NIRSpec(JWSTInstrument):

    """
    Need to over-ride get_wave_range() for NIRSpec because the effective wavelength range
    depends on the blocking filter, the aperture, and the disperser.  Also need to overload __init__
    to handle some special MSA configuration needs.
    """

    def __init__(self, mode=None, config={}, webapp=False, **kwargs):

        # Needed for 'rn' json detector parameters
        self.detector_readout_pattern = config['detector']['readout_pattern']
        config['detector']['max_total_groups'] = config['detector']['nint'] * config['detector']['ngroup']

        JWSTInstrument.__init__(self, mode=mode, config=config, webapp=webapp, **kwargs)

        slit = self.instrument.get('slit',None)
        if slit is not None:
            try:
                slit_config = self.slit_config[slit]
            except KeyError as e:
                msg = "Configuration for slit {} not specified".format(slit)
                raise DataError(value=msg)
            if self.mode == "msa":
                shutter_location = self.instrument['shutter_location']
                gap_config_file = slit_config.pop('gap')
                gap_config = read_json(os.path.join(self.ref_dir, gap_config_file), raise_except=True)
                self.shutter_locations = list(gap_config.keys())
                try:
                    self.gap = gap_config[shutter_location]
                except KeyError as e:
                    msg = "Shutter location not specified for MSA calculation: %s" % e
                    raise DataError(value=msg)
            else:
                self.gap = self.read_config_param(slit_config, "gap")

    def get_wave_range(self):
        """
        Get the wavelength range of the current instrument configuration

        Returns
        -------
        range_dict: dict
            Contains the instrument wavelength range in microns described by:
                wmin - minimum wavelength
                wmax - maximum wavelength
        """
        disperser = self.instrument['disperser']
        aperture = self.instrument['aperture']
        filt = self.instrument['filter']
        # MSA shutter configuration is read in from a separate file in a different way

        if disperser is not None:
            # get the wavelength range from the gap configuration
            if filt in self.gap[disperser]:
                # g140m and g140h have different gap configs for each blocking filter
                g_wmin = self.gap[disperser][filt]["wave_min"]
                g_wmax = self.gap[disperser][filt]["wave_max"]
            else:
                g_wmin = self.gap[disperser]["wave_min"]
                g_wmax = self.gap[disperser]["wave_max"]

            # get the wavelength range from the configuration file
            c_wmin = self.range[aperture][filt]["wmin"]
            c_wmax = self.range[aperture][filt]["wmax"]

            # get the wavelength range over which the disperser efficiency is known
            wave_blaze = self.get_wave_blaze()
            d_wmin = wave_blaze.min()
            d_wmax = wave_blaze.max()

            # get the wavelength range over which the filter throughput is known
            wave_filter = self.get_wave_filter()
            f_wmin = wave_filter.min()
            f_wmax = wave_filter.max()

            # compare filter and disperser wavelength ranges
            if f_wmax < d_wmin or d_wmax < f_wmin:
                raise RangeError(value="Disperser and Filter wavelength ranges do not overlap.")
            wmin = max(f_wmin, d_wmin, c_wmin, g_wmin)
            wmax = min(f_wmax, d_wmax, c_wmax, g_wmax)
        else:
            wmin = self.range[aperture][filt]["wmin"]
            wmax = self.range[aperture][filt]["wmax"]

        range_dict = {'wmin': wmin, 'wmax': wmax}
        return range_dict

    def create_gap_mask(self, wave):
        """
        Use the gap configuration and a wavelength vector, wave, to build a masked array that
        masks out the location of the gap.  Wavelengths between and including both gap endpoints
        will be masked out.

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to construct mask from

        Returns
        -------
        mask: numpy.ma 1D masked array
            1D array masked at the locations within the configured detector gap and 1.0 elsewhere
        """
        disperser = self.instrument['disperser']
        aperture = self.instrument['aperture']
        filt = self.instrument['filter']
        # MSA shutter configuration is read in from a separate file in a different way
        if hasattr(self, "gap"):
            gap = self.gap[disperser]

            if filt in gap:
                gap_start = gap[filt]['gap_start']
                gap_end = gap[filt]['gap_end']
            else:
                gap_start = gap['gap_start']
                gap_end = gap['gap_end']

            if gap_start is not None and gap_end is not None:
                masked_wave = ma.masked_inside(wave, gap_start, gap_end)
                mask = masked_wave / wave
                mask = ma.filled(mask, 0.0)
            else:
                mask = 1.0
        else:
            mask = 1.0
        return mask

    def get_internal_eff(self, wave):
        """
        Read in internal efficiency of NIRSpec. This is
        overloaded because the internal optical throughput is
        different for the NIRSpec IFU compared to MOS and Fixed Slit.

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        eff: numpy.ndarray or float
            Internal throughput as a function of wave
        """
        if self.mode in ["ifu", "ifu_ver"]:
            eff = self._get_throughput(wave, 'internal_ifu')
        elif self.mode in ["msa", "bots", "mos_conf", "mos_ver", "fixed_slit", "target_acq"]:
            eff = self._get_throughput(wave, 'internal_msa')
        else:
            msg = "Internal efficiency not configured for NIRSpec mode %s." % self.mode
            raise EngineInputError(value=msg)

        return eff

class NIRCam(JWSTInstrument):

    def _loadpsfs(self):
        """
        For the bar-shaped coronagraphy masks we need to load the psf_library on a per-filter basis.
        The PSF files have the aperture and filter smooshed into one string so use that as the key.
        """
        if self.instrument['aperture'] in ('masklwb', 'maskswb'):
            psf_key = "%s%s" % (self.instrument['aperture'], self.instrument['filter'])
        else:
            psf_key = self.instrument['aperture']

        self.psf_library = self._load_psf_library(psf_key)

    def get_filter_eff(self, wave):
        """
        over-ride the filter efficiency because the narrow-band filters are in the pupil wheel,
        and therefore also go through a broad-band filter in the filter wheel (which doesn't have a clear).

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        eff: numpy.ndarray or float
            Filter throughput as a function of wave
        """

        if not hasattr(self, 'double_filters'):
            msg = "NIRCam requires a mapping that describes which filters are actually a combination of two filters."
            raise DataError(value=msg)

        if self.instrument['filter'] in self.double_filters:
            eff = self._get_throughput(wave, self.instrument['filter'])
            eff_pupil = self._get_throughput(wave, self.double_filters[self.instrument['filter']])
            eff *= eff_pupil
        else:
            eff = self._get_throughput(wave, self.instrument['filter'])
        return eff

    def get_internal_eff(self, wave):
        """
        Read in internal efficiency. For NIRCam there are separate internal efficiencies for the optics common
        to all modes, the throughput of the coronagrapher substrate, and the throughputs of the optical wedges
        that bring the coronagraphy elements into the field of view of the detectors.

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        eff: numpy.ndarray or float
            Internal throughput as a function of wave
        """
        base_eff = self._get_throughput(wave, 'internal')
        coronagraphy_eff = 1.0
        wedge_eff = 1.0
        dichroic_eff = 1.0
        pupil_eff = 1.0

        # load the dichroic effects for the correct detector. Note that _get_throughput defaults to returning 1
        dichroic_eff = self._get_throughput(wave, 'dbs_{}'.format(self.instrument['detector']))

        # load the weak lens pupil transmission, if applicable
        if "wlp8" in self.instrument['aperture']:
            pupil_eff = self._get_throughput(wave, 'wlp8')

        if self.instrument['mode'] == 'coronagraphy':
            coronagraphy_eff = self._get_throughput(wave, 'coronagraphy_substrate')
            wedge_eff = self._get_throughput(wave, '{}_wedge_eff'.format(self.instrument['detector']))
        eff = base_eff * coronagraphy_eff * wedge_eff * dichroic_eff * pupil_eff

        return eff

    def dispersion_axis(self):
        """
        The dispersion axis is either along rows (the X axis) or along columns (Y axis).  By default
        it is along the X axis.  However, the GRISMC grating disperses along the Y axis as a way to
        help mitigate source crowding and confusion.

        Returns
        -------
        disp_axis: str
            Axis along with spectra are dispersed.  Allowed values are 'x' or 'y'.
        """
        if self.projection_type == 'slitless' and self.instrument['disperser'] == "grismc":
            disp_axis = "y"
        else:
            disp_axis = "x"
        return disp_axis

    def bar_width(self, x):
        """
        Width of MASKLWB or MASKSWB in arcsec as a function of X. The width at the center of the FOV is taken from the
        configuration as a function of what filter is being used.

        See NIRCam Coronagraph Operations Description, Version 4.1, Sept. 26, 2016, J.Stansberry, NIRCam Operations
        https://confluence.stsci.edu/download/attachments/52920601/NIRCam_CoronagraphOps_V4.1.pdf?version=1&modificationDate=1486495660042&api=v2
        This needs to be refactored to remove the hard-coded constants.

        Parameters
        ----------
        x: float
            X position (arcsec) in field of view

        Returns
        ------
        width: float
            Width of bar at X in arcsec
        """
        filt = self.instrument['filter']
        if (len(filt.split('_')) > 1) and (filt.split('_')[1] == 'nd'): # if this is a neutral density combination, use the main filter to get the bar width
            filt = filt.split('_')[0]
        if filt not in self.bar_offsets:
            msg = "Invalid filter, %s, for MASKLWB/MASKSWB." % filt
            raise DataError(value=msg)
        center = self.bar_offsets[filt]
        if self.instrument['aperture'] == 'maskswb':
            # the maskswb bar widens as x increases
            width = 0.2666 - 0.01777 * (center + x)
        elif self.instrument['aperture'] == 'masklwb':
            # the masklwb bar is flipped: it narrows as x increases
            center = center * -1
            x = x * -1
            width = 0.5839 - 0.03893 * (center + x)
        else:
            msg = "bar_width() method only appropriate for MASKLWB and MASKSWB apertures."
            raise EngineInputError(value=msg)
        return width

    def get_detector_qe(self, wave):
        """
        Need to over-ride get_detector_qe() to handle the two different detectors. Which one to use is keyed off
        of the configured aperture.
        """
        try:
            detector = self.instrument['detector']
        except KeyError as e:
            msg = "NIRCam aperture configuration must include which detector the aperture belongs to, sw or lw. (%s)" % e
            raise DataError(value=msg)
        qe = self._get_throughput(wave, 'qe_{}'.format(detector))

        return qe


class NIRISS(JWSTInstrument):

    """
    Need to over-ride get_internal_eff() because there are two different wheels in NIRISS.
    The 'filter wheel' has a clear position with 100% transmission. The 'pupil wheel'
    (which somewhat confusingly also contains filters) has a clear position with some
    obstruction (the Pupil Alignment Reference - PAR). So if the active filter is in the
    'filter wheel', the transmission takes a hit from the PAR in the pupil wheel.

    Also need to override dispersion_axis() because the BG150C grating disperses along the
    Y axis vs. the X axis like every other JWST disperser.
    """

    def __init__(self, mode=None, config={}, webapp=False, **kwargs):

        # these are required for calculation, but ok to live with config file defaults
        self.api_ignore = ['max_saturated_pixels', 'min_snr_threshold', 'aperture_size']

        JWSTInstrument.__init__(self, mode=mode, config=config, webapp=webapp, **kwargs)

    def dispersion_axis(self):
        """
        The dispersion axis is either along rows (the X axis) or along columns (Y axis).  By default
        it is along the X axis.  However, the GR150C grating disperses along the Y axis as a way to
        help mitigate source crowding and confusion.

        Returns
        -------
        disp_axis: str
            Axis along with spectra are dispersed.  Allowed values are 'x' or 'y'.
        """
        if self.instrument['disperser'] in ("gr150c", "gr700xd"):
            disp_axis = "y"
        else:
            disp_axis = "x"
        return disp_axis

    def get_extraction_mask(self, order):
        """
        Each SOSS order has its own extraction mask. Use the specified order to build the key to look up the
        FITS file containing the mask.

        Parameters
        ----------
        order: int
            Order whose mask to read

        Returns
        -------
        mask: 2D np.ndarray
            Mask data
        """
        if order not in (1, 2, 3):
            msg = "SOSS order %d is not valid." % order
            raise EngineInputError(value=msg)

        key = "gr700xd_%d_mask" % order

        # substrip96 is a special case where there's only one possible mask
        if self.detector['subarray'] == 'substrip96':
            key += "96"

        path = self.paths.get(key, None)

        if path is None:
            msg = "No mask configured for SOSS order %d." % order
            raise DataError(value=msg)

        mask_file = os.path.join(self.ref_dir, path)
        try:
            mask = fits.getdata(mask_file)
        except Exception as e:
            msg = "Error reading mask data for SOSS order %d: %s" % (order, type(e))
            raise DataError(value=msg)

        return mask

    def _loadpsfs(self):
        """
        Short-wavelength filters need PSFs that have the CLEAR pupil mask (0.79-2.26 microns)
        Long-wavelength filters need PSFs that have the CLEARP mask (2.37-5.04 microns)
        (See https://jwst-docs.stsci.edu/display/JTI/NIRISS+Overview, Table 2)
        Because the ranges overlap when put on a grid, we need to switch between PSF libraries.
        """
        if self.instrument['aperture'] == 'imager':
            if self.instrument['filter'] in self.lw_pupil:
                psf_key = "%s%s" % (self.instrument['aperture'], 'lw')
            else:
                psf_key = "%s%s" % (self.instrument['aperture'], 'sw')
        else:
            psf_key = self.instrument['aperture']

        self.psf_library = self._load_psf_library(psf_key)

    def get_trace(self, wave):
        """
        Read in spectral trace offsets from center of FOV. Currently wavelength-dependent spatial distortion is
        only required for SOSS mode. Other modes simply return 0's.

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate/trim trace data onto

        Returns
        -------
        trace: numpy.ndarray or float
            Spectral trace offset from center of FOV as a function of wave
        """
        if self.mode == "soss":
            disperser = self._get_disperser()
            key = "%s_wavepix" % disperser
            # handle the special case of substrip96
            if self.detector['subarray'] == 'substrip96':
                key += "96"

            # SOSS modes requires trace reference data to work properly. so raise exception if we can't load it.
            try:
                trace = self._interp_refdata(wave, key, colname='trace')
            except DataError as e:
                msg = "Spectral trace data missing for NIRISS SOSS, %s."
                raise DataError(value=msg)
        else:
            trace = np.zeros_like(wave)
        return trace

    def get_detector_pixels(self, wave):
        """
        Read in detector pixel positions for each wavelength in wave_pix

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate/trim pixel position data onto

        Returns
        -------
        pixels: numpy.ndarray or None
            Detector pixel positions corresponding to each element of wave
        """
        disperser = self._get_disperser()
        if disperser is not None:
            key = "%s_wavepix" % disperser
            # handle the special case of substrip96
            if self.detector['subarray'] == 'substrip96':
                key += "96"
            try:
                pixels = self._interp_refdata(wave, key, colname='detector_pixels', default=np.arange(len(wave)))
            except DataError as e:
                pixels = None
        else:
            pixels = None
        return pixels

class MIRI(JWSTInstrument):

    """
    The MIRI internal efficiency, detector readout patterns, etc. are more complex, and different than the other instruments,
    so some methods are redefined.
    """

    def _loadpsfs(self):
        """
        The MIRI coronagraphic target acquisition modes use a tiny box that does not include any of the obscuration, and
        only on filters that do not have the coronagraphic occulters included. Therefore, we're using different PSFs
        for them, with aperture names "fqpm1065ta", and similar.
        """
        if (self.instrument['mode'] in ('target_acq')) and (self.instrument['aperture'] in ('fqpm1065', 'fqpm1140',
                                                                                            'fqpm1550', 'lyot2300')):
            psf_key = '{0:}ta'.format(self.instrument['aperture'])
        else:
            psf_key = self.instrument['aperture']

        self.psf_library = self._load_psf_library(psf_key)

    @property
    def qe_key(self):
        """
        MIRI has three different detectors with three different QE reference files.  Use this
        method to pick the right reference file key based on the configured aperture and overload self.qe_key.

        Returns
        -------
        key: string
            Key for looking up the appropriate QE reference file
        """
        key = "miri_{}_qe".format(self.instrument['detector'])

        return key

    def get_variance_fudge(self, wave):
        """
        In addition to a scalar fudge factor, MIRI also has a chromatic variance fudge that correlates with
        the quantum yield. The information posted in Issue #2167 suggests that they want the noise scaled by an extra
        factor of the quantum yield so that the SNR scales inversely with quantum yield. The MIRI team has been asked
        to provide the chromatic fudge factor they want as a separate reference file. Until that's delivered, we'll
        use the quantum yield squared to achieve the desired effect.
        """
        scalar_fudge = self.the_detector._get_variance_fudge(wave)
        q_yield, fano_factor = self.get_quantum_yield(wave)
        var_fudge = scalar_fudge * q_yield**2

        return var_fudge

    def get_dichroic(self, type, level):
        """
        MIRI dichroic reference files are keyed by type (refl or trans), level, and disperser

        Parameters
        ----------
        type: str
            Either 'refl' for reflective or 'trans' for transmitting
        level: int
            Level 1, 2, or 3

        Returns
        -------
        key: str
            Key pointing to the dichroic efficiency data in self.paths
        """
        disperser = self.instrument['disperser']
        if disperser == 'short':
            mrs_setting = 'a'
        elif disperser == 'medium':
            mrs_setting = 'b'
        elif disperser == 'long':
            mrs_setting = 'c'

        key = "sd" + str(level) + mrs_setting + "_" + type
        return key

    def get_internal_eff(self, wave):
        """
        Calculate MIRI internal efficiency which is rather more complicated than the other instruments

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate efficiency onto

        Returns
        -------
        eff: numpy.ndarray or float
            Internal efficiency as a function of wave
        """
        aperture = self.instrument['aperture']
        mirror_eff = self.mirror_eff
        mirror_cont = self.mirror_cont
        n_refl = self.n_reflections[aperture]
        refl_eff = mirror_eff ** n_refl
        if self.mode in ['mrs', 'mrs_ts']:
            if aperture == 'ch1':
                d1r_key = self.get_dichroic('refl', 1)
                d1r_eff = self._get_throughput(wave, d1r_key)
                internal_eff = d1r_eff * refl_eff

            elif aperture == 'ch2':
                d1t_key = self.get_dichroic('trans', 1)
                d2r_key = self.get_dichroic('refl', 2)
                d1t_eff = self._get_throughput(wave, d1t_key)
                d2r_eff = self._get_throughput(wave, d2r_key)
                internal_eff = d1t_eff * d2r_eff * refl_eff

            elif aperture == 'ch3':
                d1t_key = self.get_dichroic('trans', 1)
                d2t_key = self.get_dichroic('trans', 2)
                d3r_key = self.get_dichroic('refl', 3)
                d1t_eff = self._get_throughput(wave, d1t_key)
                d2t_eff = self._get_throughput(wave, d2t_key)
                d3r_eff = self._get_throughput(wave, d3r_key)
                internal_eff = d1t_eff * d2t_eff * d3r_eff * refl_eff

            elif aperture == 'ch4':
                d1t_key = self.get_dichroic('trans', 1)
                d2t_key = self.get_dichroic('trans', 2)
                d3t_key = self.get_dichroic('trans', 3)
                d1t_eff = self._get_throughput(wave, d1t_key)
                d2t_eff = self._get_throughput(wave, d2t_key)
                d3t_eff = self._get_throughput(wave, d3t_key)
                internal_eff = d1t_eff * d2t_eff * d3t_eff * refl_eff

            else:
                raise EngineInputError(value='Invalid aperture for MIRI: {0:}'.format(aperture))

        elif (self.mode == 'coronagraphy') or (self.mode == "target_acq"):
            # The 4QPM coronagraph masks are made of germanium, and need transmission components added.
            if (aperture == 'fqpm1065') or (aperture == 'fqpm1140'):
                ar1_eff = self._get_throughput(wave, 'ge_ar1_trans')
                internal_eff = ar1_eff * refl_eff
            elif (aperture == 'fqpm1550'):
                ar2_eff = self._get_throughput(wave, 'ge_ar2_trans')
                internal_eff = ar2_eff * refl_eff
            else:
                internal_eff = np.zeros_like(wave) + refl_eff
        else:
            internal_eff = np.zeros_like(wave) + refl_eff

        # mirror contamination factor
        internal_eff = internal_eff * mirror_cont

        return internal_eff

    def get_disperser_eff(self, wave):
        """
        Overloaded here because disperser efficiency is keyed off of the aperture rather than disperser

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        disperser_eff: numpy.ndarray or float
            Disperser efficiency as a function of wave
        """
        if self.instrument['disperser'] is not None:
            key = self.instrument['aperture']
            disperser_eff = self._get_throughput(wave, key)
        else:
            disperser_eff = 1.
        return disperser_eff

    def _get_dispersion_key(self):
        """
        Overload this because the key is constructed from both the aperture and disperser rather
        than the disperser alone.

        Returns
        -------
        key: str
            Key used to get dispersion file out of self.paths
        """
        disperser = self.instrument['disperser']
        aperture = self.instrument['aperture']
        key = "%s_%s_disp" % (aperture, disperser)
        return key


def name_mapper(name=None):
    """
    General Purpose name remapping function
    If not fed a name, it returns the mapping dictionary
    If fed a name, it returns either the mapped name or the name (if no mapping is defined)
    Parameters
    ----------
    name: string
        The name of a JWST object
    Returns
    -------
    dictionary or string
        Either the remapped string (where necessary) or the complete mapping dictionary.
    """
    short_str_mappings = {
        'nircam ssgrism':  'nircam lw_tsgrism',
        'nirspec msa':     'nirspec mos',
        'ssgrism':         'lw_tsgrism',
        'msa':             'mos',
    }
    if name is None:
        return short_str_mappings
    else:
        if name in short_str_mappings:
            return short_str_mappings[name]
        else:
            return name
