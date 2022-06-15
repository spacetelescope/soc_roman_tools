# Licensed under a 3-clause BSD style license - see LICENSE.rst
# from __future__ import division, absolute_import

import copy

import numpy as np
from warnings import warn
import scipy.integrate as integrate
import scipy.interpolate as sci_int
from scipy.ndimage.interpolation import shift
import scipy.signal as sg
from astropy.io import fits

from .config import DefaultConfig
from . import astro_spectrum as astro
from . import background as bg
from . import coords
from .custom_exceptions import EngineInputError, EngineOutputError, RangeError, DataError, DataConfigurationError
from .pandeia_warnings import etc3d_warning_messages as warning_messages
from . import debug_utils




class CalculationConfig(DefaultConfig):
    """
    Encapsulate calculation configuration parameters (e.g. effects to include, noise sources to consider)
    """
    pass

class DetectorSignal(astro.ConvolvedSceneCube):
    """
    This class contains functionality for calculating the integrated electron rate of an
    astronomical source at the relevant instrument detector plane.

    Parameters
    ----------
    observation: observation.Observation instance
        Contains information required to configure a calculation
    calc_config: CalculationConfig instance
        Contains boolean flags that control which noise components are included in the calculation
    webapp: bool
        Toggle strict engine API checking
    order: int or None
        For multi-order case (i.e. SOSS), designate which order to use

    Attributes
    ----------
    """

    def __init__(self, observation, calc_config=CalculationConfig(), webapp=False, order=None, empty_scene=False):
        # Get calculation configuration
        self.calculation_config = calc_config

        # Link to the passed observation
        self.observation = observation

        # Load the instrument we're using
        self.current_instrument = observation.instrument
        # and configure it for the order we wish to use, if applicable
        self.current_instrument.order = order

        # save to the DetectorSignal instance, for convenience purposes
        self.order = order

        # how are we projecting the signal onto the detector plane?
        self.projection_type = self.current_instrument.projection_type

        # If we're in a dispersed mode, we need to know which axis the signal is dispersed along
        self.dispersion_axis = self.current_instrument.dispersion_axis()

        # Get the detector parameters (read noise, etc.)
        self.the_detector = self.current_instrument.the_detector

        # Initialize detector mask
        self.det_mask = 1.0

        # Get the background
        self.background = bg.Background(self.observation, webapp=webapp)

        # Then initialize the flux and wavelength grid
        astro.ConvolvedSceneCube.__init__(
            self,
            self.observation.scene,
            self.current_instrument,
            background=self.background,
            psf_library=self.current_instrument.psf_library,
            webapp=webapp,
            empty_scene=empty_scene
        )

        self.warnings.update(self.background.warnings)
        # We have to propagate the background through the system transmission
        # to get the background in e-/s/pixel/micron. The background rate is a 1D spectrum.
        self.bg_fp_rate = self.get_bg_fp_rate()

        # Initialize slice lists
        self.rate_list = []
        self.rate_plus_bg_list = []
        self.saturation_list = []
        self.groups_list = []
        self.types_list = []
        self.pixgrid_list = []

        # Loop over all slices and calculate the photon and electron rates through the
        # observatory for each one. Note that many modes (imaging, etc.) will have just
        # a single slice.
        for flux_cube, flux_plus_bg in zip(self.flux_cube_list, self.flux_plus_bg_list):
            # Rates for the slice without the background
            slice_rate = self.all_rates(flux_cube, add_dark_and_extended_background=False)

            # Rates for the slice with the background added
            slice_rate_plus_bg = self.all_rates(flux_plus_bg, add_dark_and_extended_background=True)

            # Saturation map for the slice
            saturation_mask = np.zeros_like(slice_rate_plus_bg['fp_pix_no_ipc_unbinned'])
            slice_saturation = self.the_detector._get_saturation_mask(slice_rate_plus_bg['fp_pix_no_ipc_unbinned'], saturation_mask)
            slice_group, sat_type = self.the_detector.get_before_sat(slice_rate_plus_bg['fp_pix_unbinned'])

            slice_saturation = self.binning.max(slice_saturation)
            slice_group = self.binning.min(slice_group)

            # The grid in the slice
            slice_pixgrid = self.get_pix_grid(slice_rate)

            # Append all slices to the master lists
            self.rate_list.append(slice_rate)
            self.rate_plus_bg_list.append(slice_rate_plus_bg)
            self.saturation_list.append(slice_saturation)
            self.groups_list.append(slice_group)
            self.types_list.append(sat_type)
            self.pixgrid_list.append(slice_pixgrid)

        # the binning operator contains the new correct spatial grid
        self.grid_unbinned = self.grid
        self.grid = self.binning.grid
        try:
            self.det_mask = self.binning.dispersion(self.det_mask)
        except TypeError:
            pass

        # Get the mapping of wavelength to pixels on the detector plane. This is grabbed from the
        # first entry in self.rate_list and is currently defined to be the same for all slices.
        self.wave_pix = self.get_wave_pix()

        # This is also grabbed from the first slice as a diagnostic
        self.fp_rate = self.get_fp_rate()

        # Reassemble rates of multiple slices on the detector
        self.rate = self.on_detector(self.rate_list)
        self.rate_unbinned = self.on_detector(self.rate_list, product_name="fp_pix_unbinned")
        self.rate_no_qyield = self.on_detector(self.rate_list, product_name="fp_pix_no_qyield")
        self.rate_plus_bg = self.on_detector(self.rate_plus_bg_list)
        self.rate_plus_bg_unbinned = self.on_detector(self.rate_plus_bg_list, product_name="fp_pix_unbinned")

        # if len(rate_list) > 1, we have a detector image multiple slices tall.
        self.binning.cropshape[0] *= len(self.rate_list)
        self.binning.newshape[0] *= len(self.rate_list)

        self.trigger_saturation_calcs()
        self.detector_pixels = self.current_instrument.get_detector_pixels(self.wave_pix)

        # JWST style
        self.brightest_pixel = np.max(self.rate_plus_bg)
        # HST Health and Safety
        self.brightest_pixel_rate = self.calc_brightest_pixel(self.rate_no_qyield,self.bg_pix_rate)
        self.detector_total_rate = self.calc_total_detector(self.rate_no_qyield,self.bg_pix_rate)

        # Get the read noise correlation matrix and store it as an attribute.
        if self.the_detector.rn_correlation:
            self.read_noise_correlation_matrix = self.current_instrument.get_readnoise_correlation_matrix(
                self.rate.shape)

    def spectral_detector_transform(self):
        """
        Create engine API format dict section containing properties of wavelength coordinates
        at the detector plane.

        Returns
        -------
        t: dict (engine API compliant keys)
        """
        t = {}
        t['wave_det_refpix'] = 0
        t['wave_det_max'] = self.wave_pix.max()
        t['wave_det_min'] = self.wave_pix.min()

        # there are currently three projection_type's which are basically detector plane types:
        #
        # 'spec' - where the detector plane is purely dispersion vs. spatial
        # 'slitless' - basically a special case of 'spec' with where dispersion and spatial are mixed
        # 'image' - where the detector plane is purely spatial vs. spatial (i.e. no disperser element)
        #
        # 'IFU' mode is of projection_type='spec' because the mapping from detector X pixels to
        # wavelength is the same for each slice.  this projection_type will work for 'MSA' mode as well
        # because we will only handle one aperture at a time.  'slitless' spectroscopy will mix
        # spatial and dispersion information onto the detector X axis.  however, the detector
        # plane is fundamentally spatial vs. wavelength in that case so it's handled the same as
        # projection_type='spec'. creating a spectrum for a specific target will be handled via the
        # extraction strategy.
        if self.projection_type in ('spec', 'slitless', 'multiorder'):
            t['wave_det_size'] = len(self.wave_pix)
            if len(self.wave_pix) > 1:
                # we don't yet have a way of handling non-linear coordinate transforms here. that said,
                # this is mostly right for most of our cases with nirspec prism being the notable exception.
                # this is also only used for plotting purposes while the true actual wave_pix mapping is used
                # internally for all calculations.
                t['wave_det_step'] = (self.wave_pix[-1] - self.wave_pix[0]) / t['wave_det_size']
            else:
                t['wave_det_step'] = 0.0
            t['wave_det_refval'] = self.wave_pix[0]
        elif self.projection_type == "image":
            t['wave_det_step'] = 0.0
            t['wave_det_refval'] = self.wave_pix[0]
            t['wave_det_size'] = 1
        else:
            message = "Unsupported projection_type: %s" % self.projection_type
            raise EngineOutputError(value=message)
        return t

    def wcs_info(self):
        """
        Get detector coordinate transform as a dict of WCS keyword/value pairs.

        Returns
        -------
        header: dict
            WCS header keys defining coordinate transform in the detector plane
        """
        if self.projection_type == 'image':
            # if we're in imaging mode, the detector sampling is the same as the model
            header = self.grid.wcs_info()
        elif self.projection_type in ('spec', 'slitless', 'multiorder'):
            # if we're in a dispersed mode, dispersion can be either along the X or Y axis. the image outputs in
            # the engine Report are rotated so that dispersion will always appear to be along the X axis with
            # wavelength increasing with increasing X (i. e. dispersion angle of 0).  currently, the only other
            # supported dispersion angle is 90 which is what we get when dispersion_axis == 'y'.
            t = self.grid.as_dict()
            t.update(self.spectral_detector_transform())
            header = {
                'ctype1': 'Wavelength',
                'crpix1': 1,
                'crval1': t['wave_det_min'] - 0.5 * t['wave_det_step'],
                'cdelt1': t['wave_det_step'],
                'cunit1': 'um',
                'cname1': 'Wavelength',
                'ctype2': 'Y offset',
                'crpix2': 1,
                'crval2': t['y_min'] - 0.5 * t['y_step'],
                'cdelt2': -t['y_step'],
                'cunit2': 'arcsec',
                'cname2': 'Detector Offset',
            }
            if self.dispersion_axis == 'y':
                header['ctype2'] = 'X offset'
                header['crval2'] = t['x_min'] - 0.5 * t['x_step'],
                header['cdelt2'] = t['x_step']
        else:
            message = "Unsupported projection_type: %s" % self.projection_type
            raise EngineOutputError(value=message)
        return header

    def get_wave_pix(self):
        """
        Return the mapping of wavelengths to pixels on the detector plane
        """
        return self.rate_list[0]['wave_pix']

    def get_fp_rate(self):
        """
        Return scene flux at the focal plane in e-/s/pixel/micron (excludes background)
        """
        return self.rate_list[0]['fp']

    def get_bg_fp_rate(self):
        """
        Calculate background in e-/s/pixel/micron at the focal plane. Also correct for any excess in predicted background
        if there are pupil losses in the PSF. (#2529)
        """
        bg_fp_rate = self.focal_plane_rate(self.ote_rate(self.background.mjy_pix))
        wave_range = self.current_instrument.get_wave_range()
        pupil_thru = self.current_instrument.psf_library.get_pupil_throughput(wave_range['wmin'],
                                                                              self.current_instrument.instrument[
                                                                                  'instrument'],
                                                                              self.current_instrument.instrument[
                                                                                  'aperture'])
        return bg_fp_rate * pupil_thru

    def get_bg_pix_rate(self):
        """
        Calculate the background on the detector in e-/s/pixel

        By computing off of these on-detector stitched properties, we can check if ANY of the background is saturating.
        """
        bg_pix_rate = self.rate_plus_bg_unbinned - self.rate_unbinned
        return bg_pix_rate

    def on_detector(self, rate_list, product_name="fp_pix"):
        """
        This will take the list of (pixel) rates and use them create a single detector frame. A single
        image will only have one rate in the list, but the IFUs will have n_slices. There may be other examples,
        such as different spectral orders for NIRISS. It is not yet clear how many different flavors there are, so
        this step may get refactored if it gets too complicated. Observing modes that only have one set of rates
        (imaging and single-slit spectroscopy, for instance) will still go through this, but the operation is trivial.
        """
        aperture_sh = rate_list[0][product_name].shape
        n_apertures = len(rate_list)
        detector_shape = (aperture_sh[0] * n_apertures, aperture_sh[1])
        detector = np.zeros(detector_shape)

        i = 0
        for rate in rate_list:
            detector[i * aperture_sh[0]:(i + 1) * aperture_sh[0], :] = rate[product_name]
            i += 1

        return detector

    def get_pix_grid(self, rate):
        """
        Generate the coordinate grid of the detector plane
        """
        if self.projection_type == 'image':
            grid = self.grid
        elif self.projection_type in ('spec', 'slitless', 'multiorder'):
            nw = rate['wave_pix'].shape[0]
            if self.dispersion_axis == 'x':
                # for slitless calculations, the dispersion axis is longer than the spectrum being dispersed
                # because the whole field of view is being dispersed. 'excess' is the size of the FOV
                # and half will be to the left of the blue end of the spectrum and half to the right of the red end.
                # this is used to create the new spatial coordinate transform for the pixel image on the detector.
                excess = rate['fp_pix'].shape[1] - nw
                pix_grid = coords.IrregularGrid(
                    self.grid.col,
                    (np.arange(nw + excess) - (nw + excess) / 2.0) * self.grid.xsamp
                )
            else:
                excess = rate['fp_pix'].shape[0] - nw
                pix_grid = coords.IrregularGrid(
                    (np.arange(nw + excess) - (nw + excess) / 2.0) * self.grid.ysamp,
                    self.grid.row
                )
            return pix_grid
        else:
            raise EngineOutputError(value="Unsupported projection_type: %s" % self.projection_type)
        return grid

    def all_rates(self, flux, add_dark_and_extended_background=False):
        """
        Calculate rates in e-/s/pixel/micron or e-/s/pixel given a flux cube in mJy

        Parameters
        ----------
        flux: ConvolvedSceneCube instance
            Convolved source flux cube with flux units in mJy
        add_dark_and_extended_background: bool (default=False)
            Toggle for including extended background not contained within the flux cube, and
            also the dark current

        Returns
        -------
        products: dict
            Dict of products produced by rate calculation.
                'wave_pix' - Mapping of wavelength to detector pixels
                'ote' - Source rate at the telescope aperture
                'fp' - Source rate at the focal plane in e-/s/pixel/micron
                'fp_pix' - Source rate per pixel
                'fp_pix_no_ipc' - Source rate per pixel excluding effects if inter-pixel capacitance
        """
        # The source rate at the telescope aperture
        ote_rate = self.ote_rate(flux)

        # The source rate at the focal plane in interacting photons/s/pixel/micron
        fp_rate = self.focal_plane_rate(ote_rate)

        total_rate = fp_rate.sum(axis=(0,1))

        # the fp_pix_variance is the variance of the per-pixel electron rate and includes the chromatic effects
        # of quantum yield.
        if self.projection_type == 'image':
            # The wavelength-integrated rate in e-/s/pixel, relevant for imagers
            fp_pix_rate, fp_pix_no_qyield, fp_pix_variance = self.image_rate(fp_rate)
            wave_pix = self.wave_eff(fp_rate)

        elif self.projection_type == 'spec':
            # The wavelength-integrated rate in e-/s/pixel, relevant for spectroscopy
            wave_pix, fp_pix_rate, fp_pix_no_qyield, fp_pix_variance = self.spec_rate(fp_rate)

        elif self.projection_type in ('slitless', 'multiorder'):
            # The wavelength-integrated rate in e-/s/pixel, relevant for slitless spectroscopy
            wave_pix, fp_pix_rate, fp_pix_no_qyield, fp_pix_variance = self.slitless_rate(
                fp_rate,
                add_extended_background=add_dark_and_extended_background
            )

        else:
            raise EngineOutputError(value="Unsupported projection_type: %s" % self.projection_type)

        # Dark is added if and only if we're processing a rate with a background
        if add_dark_and_extended_background:
            fp_pix_rate += self.the_detector.dark_current
            fp_pix_variance += self.the_detector.dark_current

        # Include IPC effects, if available and requested
        if self.the_detector.ipc:
            kernel = self.current_instrument.get_ipc_kernel()
            fp_pix_rate_ipc = self.ipc_convolve(fp_pix_rate, kernel)
        else:
            fp_pix_rate_ipc = copy.deepcopy(fp_pix_rate)

        fp_pix_unbinned = fp_pix_rate_ipc
        fp_pix_no_ipc_unbinned = fp_pix_rate
        fp_pix_no_qyield_unbinned = fp_pix_no_qyield
        fp_pix_variance_unbinned = fp_pix_variance
        wave_pix = self.binning.dispersion(wave_pix)
        fp_pix_rate_ipc = self.binning.sum(fp_pix_rate_ipc)
        fp_pix_rate = self.binning.sum(fp_pix_rate)
        fp_pix_no_qyield = self.binning.sum(fp_pix_no_qyield)
        fp_pix_variance = self.binning.sum(fp_pix_variance)

        # fp_pix is the final product. Since there is no reason to
        # carry around the ipc label everywhere, we rename it here.
        products = {
            'wave_pix': wave_pix,
            'ote': ote_rate,
            'fp': fp_rate,
            'fp_pix': fp_pix_rate_ipc,
            'fp_pix_no_ipc': fp_pix_rate,  # this is for calculating saturation
            'fp_pix_no_qyield': fp_pix_no_qyield, # this is for calculating feasibility thread
            'fp_pix_variance': fp_pix_variance,  # this is for calculating the detector noise
            'fp_pix_unbinned': fp_pix_unbinned, # this is fp_pix before binning
            'fp_pix_no_ipc_unbinned': fp_pix_no_ipc_unbinned, # this is fp_pix_no_ipc before binning
            'fp_pix_no_qyield_unbinned': fp_pix_no_qyield_unbinned, # this is calculating feasibility thread 
            'fp_pix_variance_unbinned': fp_pix_variance_unbinned  # this is for calculating the detector noise
        }
        return products

    def ote_rate(self, flux):
        """
        Calculate source rate in e-/s/pixel/micron at the telescope entrance aperture given
        a flux cube in mJy/pixel.
        """
        # spectrum in mJy/pixel, wave in micron, f_lambda in photons/cm^2/s/micron
        f_lambda = 1.5091905 * (flux / self.wave)
        ote_int = self.current_instrument.telescope.get_ote_eff(self.wave)
        coll_area = self.current_instrument.telescope.coll_area
        a_lambda = coll_area * ote_int
        # e-/s/pixel/micron
        ote_rate = f_lambda * a_lambda
        return ote_rate

    def focal_plane_rate(self, rate):
        """
        Takes the output from self.ote_rate() and multiplies it by the components
        of efficiency within the system and returns the source rate at the focal plane in
        e-/s/pixel/micron.
        """
        filter_eff = self.current_instrument.get_filter_eff(self.wave)
        disperser_eff = self.current_instrument.get_disperser_eff(self.wave)
        internal_eff = self.current_instrument.get_internal_eff(self.wave)
        qe = self.current_instrument.get_detector_qe(self.wave)

        fp_rate = rate * filter_eff * disperser_eff * internal_eff * qe
        return fp_rate

    def spec_rate(self, rate):
        '''
        For slitted spectrographs, calculate the detector signal by integrating
        along the dispersion direction of the cube (which is masked by a, by assumption,
        narrow slit). For slitless systems or slits wider than the PSF, the slitless_rate
        method should be used to preserve spatial information within the slit.

        Parameters
        ---------
        rate: numpy.ndarray
            Rate of photons interacting with detector as a function of model wavelength set

        Returns
        -------
        products: 3-element tuple of numpy.ndarrays
            first element - map of pixel to wavelength
            second element - electron rate per pixel
            third element - variance of electron rate per pixel
        '''
        dispersion = self.current_instrument.get_dispersion(self.wave)
        wave_pix = self.current_instrument.get_wave_pix()
        wave_pix_trunc = wave_pix[np.where(np.logical_and(wave_pix >= self.wave.min(),
                                                          wave_pix <= self.wave.max()))]

        # Check that the source spectrum is actually inside the instrumental wavelength
        # coverage.
        if len(wave_pix_trunc) == 0:
            raise RangeError(value='wave and wave_pix do not overlap')

        # Check the dispersion axis to determine which axis to sum and interpolate over
        if self.dispersion_axis == 'x':
            axis = 1
        else:
            axis = 0

        # We can simply sum over the dispersion direction. This is where we lose the spatial information within the aperture.
        spec_rate = np.sum(rate, axis=axis)

        # And then scale to the dispersion function (pixel/micron) to transform
        # from e-/s/micron to e-/s/pixel.
        spec_rate_pix = spec_rate * dispersion

        # but we are still sampled on the internal grid, so we have to interpolate to the pixel grid.
        # use kind='slinear' since it's ~2x more memory efficient than 'linear'. 'slinear' uses different code path to
        # calculate the slopes.
        int_spec_rate = sci_int.interp1d(self.wave, spec_rate_pix, axis=axis, kind='slinear', assume_sorted=True,
                                         copy=False)
        spec_rate_pix_sampled = int_spec_rate(wave_pix_trunc)

        # Handle a detector gap here by constructing a mask. If the current_instrument implements it,
        # it'll be a real mask array.  Otherwise it will simply be 1.0.
        self.det_mask = self.current_instrument.create_gap_mask(wave_pix_trunc)

        # this is the interacting photon rate in the detector with mask applied.
        spec_rate_pix_sampled *= self.det_mask

        # Add effects of non-unity quantum yields. For the spec projection, we assume that the quantum yield does not
        # change over a spectral element. Then we can just multiply the products by the relevant factors.
        q_yield, fano_factor = self.current_instrument.get_quantum_yield(wave_pix_trunc)

        # convert the photon rate to electron rate by multiplying by the quantum yield which is a function of wavelength
        spec_electron_rate_pix = spec_rate_pix_sampled * q_yield
        spec_rate_no_qyield = spec_rate_pix_sampled

        # to meet IDT expectations, some instruments require a possibly chromatic fudge factor to be applied
        # to the per-pixel electron rate variance.
        var_fudge = self.current_instrument.get_variance_fudge(wave_pix_trunc)

        # the variance in the electron rate, Ve, is also scaled by the quantum yield plus a fano factor which is
        # analytic in the simple 1 or 2 electron case: Ve = (qy + fano) * Re.  since Re is the photon rate
        # scaled by the quantum yield, Re = qy * Rp, we get: Ve = qy * (qy + fano) * Rp
        spec_electron_variance_pix = spec_rate_pix_sampled * q_yield * (q_yield + fano_factor) * var_fudge

        # bin the instrument wavelength grid in pixels
        if self.dispersion_axis == 'x':
            bin_x = self.current_instrument.get_dispersion_binning()
            bin_y = self.current_instrument.get_spatial_binning()
        else:
            bin_x = self.current_instrument.get_spatial_binning()
            bin_y = self.current_instrument.get_dispersion_binning()
        self.binning = DetectorBinning(spec_electron_rate_pix, bin_x, bin_y, self.grid, dispersion_axis=self.dispersion_axis, wave=wave_pix_trunc)

        products = wave_pix_trunc, spec_electron_rate_pix, spec_rate_no_qyield, spec_electron_variance_pix

        return products

    def image_rate(self, rate):
        '''
        Calculate the electron rate for imaging modes by integrating along
        the wavelength direction of the cube.

        Parameters
        ---------
        rate: numpy.ndarray
            Rate of photons interacting with detector as a function of model wavelength set

        Returns
        -------
        products: 2-element tuple of numpy.ndarrays
            first element - electron rate per pixel
            second element - variance of electron rate per pixel
        '''
        q_yield, fano_factor = self.current_instrument.get_quantum_yield(self.wave)

        # convert the photon rate to electron rate by multiplying by the quantum yield which is a function of wavelength
        electron_rate_pix = integrate.simps(rate * q_yield, self.wave)
        image_rate_no_qyield = integrate.simps(rate, self.wave)

        # to meet IDT expectations, some instruments require a possibly chromatic fudge factor to be applied
        # to the per-pixel electron rate variance.
        var_fudge = self.current_instrument.get_variance_fudge(self.wave)

        # the variance in the electron rate, Ve, is also scaled by the quantum yield plus a fano factor which is
        # analytic in the simple 1 or 2 electron case: Ve = (qy + fano) * Re.  since Re is the photon rate
        # scaled by the quantum yield, Re = qy * Rp, we get: Ve = qy * (qy + fano) * Rp
        electron_variance_pix = integrate.simps(rate * q_yield * (q_yield + fano_factor) * var_fudge, self.wave)

        # bin the instrument wavelength grid in pixels
        bin_x = self.current_instrument.get_spatial_binning()
        bin_y = self.current_instrument.get_spatial_binning()
        self.binning = DetectorBinning(electron_rate_pix, bin_x, bin_y, self.grid)

        products = electron_rate_pix, image_rate_no_qyield, electron_variance_pix

        return products

    def slitless_rate(self, rate, add_extended_background=True):
        '''
        Calculate the detector rates for slitless modes. Here we retain all spatial information and build
        up the detector plane by shifting and coadding the frames from the convolved flux cube. Also need to handle
        and add background that comes from outside the flux cube, but needs to be accounted for.

        How this works:
        1. We calculate wave_pix and truncate it to the wavelength range of the instrument (which should already 
        have been done, but it is a check)
        2. We grab the detector_pixel (dispersion) and trace (cross-dispersion) pixel locations, interpolated per
        wavelength on the detector. Except for modes with specific tracefile that map wavelength to 
        (detector_pixel,trace), these will simply be numerical  and constant-0 files specifying the upper left 
        corner of where to put the wavelength grid.
        3. We then step through wavelengths. If there IS a trace file, we shift in the cross dispersion direction
        and add the rate cube plane to the image at that location. Otherwise, no expensive shifting operation is
        necessary. 
        4. If extended background is switched on, we add background for that wavelength to the rest of the detector
        image not covered by the rate cube plane, because background from all wavelengths will fall on all parts of
        a slitless calculation.

        Parameters
        ----------
        rate: 3D numpy.ndarray
            Cube containing the flux rate at the focal plane
        add_extended_background: bool (default: True)
            Toggle for including extended background not contained within the flux cube

        Returns
        -------
        products: 2 entry tuple
            wave_pix: 1D numpy.ndarray containing wavelength to pixel mapping on the detector plane
            spec_rate: 2D numpy.ndarray of detector count rates
        '''
        wave_pix = self.current_instrument.get_wave_pix()
        wave_subs = np.where(
            np.logical_and(
                wave_pix >= self.wave.min(),
                wave_pix <= self.wave.max()
            )
        )
        wave_pix_trunc = wave_pix[wave_subs]

        if len(wave_pix_trunc) == 0:
            raise RangeError(value='wave and wave_pix do not overlap')

        dispersion = self.current_instrument.get_dispersion(wave_pix_trunc)
        trace = self.current_instrument.get_trace(wave_pix_trunc)
        dpix = self.current_instrument.get_detector_pixels(wave_pix_trunc).astype(np.int)

        if dispersion.shape != wave_pix_trunc.shape:
            raise DataConfigurationError("Dispersion and wavelength axes do not match! This should never happen.")

        q_yield, fano_factor = self.current_instrument.get_quantum_yield(wave_pix_trunc)
        # if we kind='slinear' since it's ~2x more memory efficient than 'linear'. 'slinear' uses different code
        # path to calculate the slopes. However, slinear is *much* slower, so it is a tradeoff. Also lowering the
        # rate type to float32 to conserve memory.
        int_rate_pix = sci_int.interp1d(self.wave, rate.astype(np.float32, casting='same_kind'),
                                        kind='linear', axis=2, assume_sorted=True, copy=False)
        rate_pix = int_rate_pix(wave_pix_trunc)

        # convert the photon rate to electron rate by multiplying by the quantum yield which is a function of wavelength
        electron_rate_pix = rate_pix * q_yield
        electron_rate_no_qyield = rate_pix

        # to meet IDT expectations, some instruments require a possibly chromatic fudge factor to be applied
        # to the per-pixel electron rate variance.
        var_fudge = self.current_instrument.get_variance_fudge(wave_pix_trunc)

        # the variance in the electron rate, Ve, is also scaled by the quantum yield plus a fano factor which is
        # analytic in the simple 1 or 2 electron case: Ve = (qy + fano) * Re.  since Re is the photon rate
        # scaled by the quantum yield, Re = qy * Rp, we get: Ve = qy * (qy + fano) * Rp
        electron_variance_pix = rate_pix * q_yield * (q_yield + fano_factor) * var_fudge

        # interpolate the background onto the pixel spacing
        int_bg_fp_rate = sci_int.interp1d(self.wave, self.bg_fp_rate.astype(np.float32, casting='same_kind'),
                                          kind='linear', assume_sorted=True, copy=False)
        bg_fp_rate_pix = int_bg_fp_rate(wave_pix_trunc)

        # calculate electron rate and variance due to background
        bg_electron_rate = bg_fp_rate_pix * q_yield
        bg_electron_variance = bg_fp_rate_pix * q_yield * (q_yield + fano_factor) * var_fudge

        # The first part of this code is meant to add the PSF images from the convolved scene cube along either the x
        # or y axis depending on the dispersion axis, optionally following the path of a spectral trace (currently used
        # only for SOSS mode). The psfs will be added to all locations within the resolution element.
        #
        # Because, in slitless modes, the disperser is dispersing light coming in from everywhere in the pupil plane,
        # every part of the detector should have a contribution from every wavelength of light (from both orders, for
        # SOSS mode). The add_extended_background statement does that - it fills every pixel up to i, and after
        # i+rate_pix.shape[1], with the same background that comes baked into the rate_pix images thanks to the
        # AdvancedPSF functions that create the convolved scene cube.

        if self.empty_scene:
            # if we have an explicitly empty scene, we're doing a background-only order calculation and don't need to
            # even pretend to disperse the spectrum - CombinedSignal will handle padding it to match the interesting
            # order(s) and this way there will be no need to trim.
            spec_shape = (rate_pix.shape[0], rate_pix.shape[1])
            spec_rate = np.zeros(spec_shape)
            spec_rate_no_qyield = np.zeros(spec_shape)
            spec_variance = np.zeros(spec_shape)
            if add_extended_background:
                for i in np.arange(dispersion.shape[0]):
                    spec_rate += bg_electron_rate[i] * dispersion[i]
                    spec_rate_no_qyield += bg_fp_rate_pix[i] * dispersion[i]
                    spec_variance += bg_electron_variance[i] * dispersion[i]
        else:  # if the scene is data
            # dispersion_axis tells us whether we need to sum the planes of the cube horizontally
            # or vertically on the detector plane.
            if self.dispersion_axis == 'x':
                spec_shape = (rate_pix.shape[0], rate_pix.shape[2] + rate_pix.shape[1] + dpix[0])
                spec_rate = np.zeros(spec_shape)
                spec_rate_no_qyield = np.zeros(spec_shape)
                spec_variance = np.zeros(spec_shape)
                for i in range(len(wave_pix_trunc)):
                    # Background not yet completely added. Make sure there is a trace shift to be done so that we
                    # don't make an expensive call to shift() if we don't have to. Use mode='nearest' to fill in new
                    # pixels with background when image is shifted.
                    # dispersion pixel location for any given wavelength is dpix, cross-dispersion pixel location is trace
                    if trace[i] != 0.0:
                        spec_rate[:, dpix[i]:dpix[i] + rate_pix.shape[1]] += shift(
                            electron_rate_pix[:, :, i],
                            shift=(trace[i], 0),
                            mode='nearest',
                            order=1
                        ) * dispersion[i]
                        spec_rate_no_qyield[:, dpix[i]:dpix[i] + rate_pix.shape[1]] += shift(
                            electron_rate_no_qyield[:, :, i],
                            shift=(trace[i], 0),
                            mode='nearest',
                            order=1
                        ) * dispersion[i]
                        spec_variance[:, dpix[i]:dpix[i] + rate_pix.shape[1]] += shift(
                            electron_variance_pix[:, :, i],
                            shift=(trace[i], 0),
                            mode='nearest',
                            order=1
                        ) * dispersion[i]
                    else:
                        spec_rate[:, dpix[i]:dpix[i] + rate_pix.shape[1]] += electron_rate_pix[:, :, i] * dispersion[i]
                        spec_rate_no_qyield[:, dpix[i]:dpix[i] + rate_pix.shape[1]] += electron_rate_no_qyield[:, :, i] * dispersion[i]
                        spec_variance[:, dpix[i]:dpix[i] + rate_pix.shape[1]] += electron_variance_pix[:, :, i] * dispersion[i]

                    # Adding background to all other pixels, unless we are asked not to.
                    if add_extended_background:
                        spec_rate[:, :dpix[i]] += bg_electron_rate[i] * dispersion[i]
                        spec_rate[:, dpix[i] + rate_pix.shape[1]:] += bg_electron_rate[i] * dispersion[i]
                        spec_rate_no_qyield[:, :dpix[i]] += bg_fp_rate_pix[i] * dispersion[i]
                        spec_rate_no_qyield[:, dpix[i] + rate_pix.shape[1]:] += bg_fp_rate_pix[i] * dispersion[i]
                        spec_variance[:, :dpix[i]] += bg_electron_variance[i] * dispersion[i]
                        spec_variance[:, dpix[i] + rate_pix.shape[1]:] += bg_electron_variance[i] * dispersion[i]
            else:  # if the dispersion is on the y axis
                spec_shape = (rate_pix.shape[2] + rate_pix.shape[0] + dpix[0], rate_pix.shape[1])
                spec_rate = np.zeros(spec_shape)
                spec_rate_no_qyield = np.zeros(spec_shape)
                spec_variance = np.zeros(spec_shape)
                for i in range(len(wave_pix_trunc)):
                    # Background not yet completely added. Make sure there is a trace shift to be done so that we
                    # don't make an expensive call to shift() if we don't have to. Use mode='nearest' to fill in new
                    # pixels with background when image is shifted.
                    # dispersion pixel location for any given wavelength is dpix, cross-dispersion pixel location is trace
                    if trace[i] != 0.0:
                        spec_rate[dpix[i]:dpix[i] + rate_pix.shape[0], :] += shift(
                            electron_rate_pix[:, :, i],
                            shift=(0, trace[i]),
                            mode='nearest',
                            order=1
                        ) * dispersion[i]
                        spec_rate_no_qyield[dpix[i]:dpix[i] + rate_pix.shape[0], :] += shift(
                            electron_rate_no_qyield[:, :, i],
                            shift=(0, trace[i]),
                            mode='nearest',
                            order=1
                        ) * dispersion[i]
                        spec_variance[dpix[i]:dpix[i] + rate_pix.shape[0], :] += shift(
                            electron_variance_pix[:, :, i],
                            shift=(0, trace[i]),
                            mode='nearest',
                            order=1
                        ) * dispersion[i]
                    else:
                        spec_rate[dpix[i]:dpix[i] + rate_pix.shape[0], :] += electron_rate_pix[:, :, i] * dispersion[i]
                        spec_rate_no_qyield[dpix[i]:dpix[i] + rate_pix.shape[0], :] += electron_rate_no_qyield[:, :, i] * dispersion[i]
                        spec_variance[dpix[i]:dpix[i] + rate_pix.shape[0], :] += electron_variance_pix[:, :, i] * dispersion[i]
                    # Adding background to all other pixels, unless we are asked not to.
                    if add_extended_background:
                        spec_rate[:dpix[i], :] += bg_electron_rate[i] * dispersion[i]
                        spec_rate[dpix[i] + rate_pix.shape[0]:, :] += bg_electron_rate[i] * dispersion[i]
                        spec_rate_no_qyield[:dpix[i], :] += bg_fp_rate_pix[i] * dispersion[i]
                        spec_rate_no_qyield[dpix[i] + rate_pix.shape[0]:, :] += bg_fp_rate_pix[i] * dispersion[i]
                        spec_variance[:dpix[i], :] += bg_electron_variance[i] * dispersion[i]
                        spec_variance[dpix[i] + rate_pix.shape[0]:, :] += bg_electron_variance[i] * dispersion[i]

        # bin the instrument wavelength grid in pixels
        if self.dispersion_axis == 'x':
            bin_x = self.current_instrument.get_dispersion_binning()
            bin_y = self.current_instrument.get_spatial_binning()
        else:
            bin_x = self.current_instrument.get_spatial_binning()
            bin_y = self.current_instrument.get_dispersion_binning()
        self.binning = DetectorBinning(spec_rate, bin_x, bin_y, self.grid, dispersion_axis=self.dispersion_axis, wave=wave_pix_trunc, projection_type=self.projection_type)

        # dispersion_axis determines whether wavelength is the first or second axis
        if self.dispersion_axis == 'x' or self.projection_type == 'multiorder':
            products = wave_pix_trunc, spec_rate, spec_rate_no_qyield, spec_variance
        else:
            # if dispersion is along Y, wavelength increases bottom to top, but Y index increases top to bottom.
            # flip the Y axis to account for this.
            products = wave_pix_trunc, np.flipud(spec_rate), np.flipud(spec_rate_no_qyield), np.flipud(spec_variance)

        return products

    def wave_eff(self, rate):
        rate_tot = np.nansum(rate, axis=0)
        a = np.sum(rate_tot * self.wave)
        b = np.sum(rate_tot)
        if (b > 0.0) and (a > 0.0):
            wave_eff = a / b
        else:
            wave_eff = self.wave.mean()
        wave_eff_arr = np.array([wave_eff])
        return wave_eff_arr

    def get_projection_type(self):
        return self.projection_type

    def ipc_convolve(self, rate, kernel):
        fp_pix_ipc = sg.fftconvolve(rate, kernel, mode='same')

        debug_utils.debugarrays.store('signal', 'ipc_convolve',
                                      {
                                          'rate': rate,
                                          'kernel': kernel,
                                          'fp_pix_ipc': fp_pix_ipc,
                                          'description': 'This is just a short, unnecessary description.'
                                      })

        return fp_pix_ipc

    def get_saturation_mask(self, rate=None):
        """
        Return a numpy array indicating pixels with full saturation (2), partial saturation (1) and no saturation (0).

        Parameters
        ----------
        rate: None or 2D np.ndarray
            Detector plane rate image used to build saturation map from

        Returns
        -------
        mask: 2D np.ndarray
            Saturation mask image
        """
        if rate is None:
            rate = self.rate_plus_bg_unbinned

        saturation_mask = np.zeros(rate.shape)

        # if (a if condition else b)
        # if the calculation_config value is True or False (edited by user), use True or False. If it's none, use the instrument-team-defined default from detector_config
        if (self.calculation_config.effects['saturation'] if self.calculation_config.effects['saturation'] is not None else self.the_detector.saturation):
            saturation_mask = self.the_detector._get_saturation_mask(rate, saturation_mask)

        saturation_mask = self.binning.max(saturation_mask)

        return saturation_mask

    def calc_total_detector(self, rate, bg_pix_rate):
        """ 
        Compute the total detector rate.

        Pandeia only handles a small "postage stamp" FOV, which means any given
        DetectorSignal only covers a portion of the entire detector.

        The total detector rate is built out of four rates:
        1. The source signal (from rate_no_qyield). Should be entirely 
           contained in that computed 2D projection.
        2. The sky (and thermal) background flux. If there's a slit, it should be
           contained entirely within the slit height and dispersed width, and it has 
           wavelength dependence in the dispersion direction. If there's no slit, it 
           covers the entire exposed portion of the detector and is gray.
        3. The dark current. Covers the entire detector, exposed or not.
        4. Scattered light, if defined. Covers the same region as the sky background.

        Other assumptions:
        The 2D projection may be a focal plane consisting of more than one detector. 
        In this case, the dark_current is not a scalar, but a 2D array giving the
        spatially variable dark current that will have to be cropped to the single
        detector.

        The background we have is a background that CONTAINS the dark current, which 
        must be subtracted off before further processing. We then have a clean background
        over the area of the "postage stamp" which we can scale accordingly. In the case 
        of a slit, even though the background is not gray it's only being extended in one
        direction, so we can still simply multiply the background flux by a scaling factor.
        The outlying case is a slit whose height is less than the postage stamp, in which
        case it's fully contained by the bg_pix_rate array and no scaling is needed.

        """

        chips = self.current_instrument.get_chip_dimensions()

        detector_total_rate_dict = {}
        for chip in chips:
            ## Process the source signal
            rate_temp = rate[chip["bounds"]["y"],chip["bounds"]["x"]]
            rate_full = np.sum(rate_temp)

            ## Process the dark current
            # dark current is either a scalar (float) or 2D map
            if isinstance(self.the_detector.dark_current, float):
                dc = self.the_detector.dark_current
                dc_full = chip["size"]["y"] * chip["size"]["x"] * dc
            else:
                dc = self.the_detector.dark_current[chip["bounds"]["y"],chip["bounds"]["x"]]
                dc_full = np.sum(dc)

            ## Process the background
            bg = bg_pix_rate[chip["bounds"]["y"],chip["bounds"]["x"]] 
            bg -= np.mean(dc) # dc is either a scalar or the entire size of the detector
            if self.projection_type == "spec":
                if self.dispersion_axis == "x":
                    bgx = rate_temp.shape[1]
                    # if this is fractional, it's actually fine. All we need is the number, anyway.
                    bgy = self.current_instrument.get_slit_pars()['xdisp'] / \
                               self.current_instrument.get_aperture_pars()["plate_scale"][0]
                    
                else:
                    bgx = self.current_instrument.get_slit_pars()['xdisp'] / \
                               self.current_instrument.get_aperture_pars()["plate_scale"][1]
                    bgy = rate_temp.shape[0]
                
            else:
                bgx = chip["size"]["x"]
                bgy = chip["size"]["y"]
            # if either slit-limited dimension is SMALLER than the background array, the background
            # array contains the entire background flux and needs no scaling.
            if bgx < chip["size"]["x"] or bgy < chip["size"]["y"]:
                bg_full = np.sum(bg)
            else:
                # we are multiplying the background by a factor of how much background is covered by the slit
                bg_full = np.mean(bg) * bgy * bgx

            total = rate_full + dc_full + bg_full 

            # selection logic: If calculation_config specifies it, do that. If not, use the predefined data value
            scatter = scatter = self.calculation_config.noise['scatter'] if self.calculation_config.noise['scatter'] is not None else self.the_detector.scatter
            if scatter:
                ## Process scattered light
                # Scattering depends on the disperser and cenwave
                global_scattering = self.current_instrument.get_global_scattering()
                scattered_full = global_scattering * rate_full

                total += scattered_full

            detector_total_rate_dict[chip["name"]] = total

        return detector_total_rate_dict

    def calc_brightest_pixel(self, rate, bg_pix_rate):
        """
        Compute the brightest pixel rate.

        Pandeia only handles a small "postage stamp" FOV, which means any given
        DetectorSignal only covers a portion of the entire detector. Using the 
        assumption that all the source flux is contained within that postage stamp,
        we can extend the background and dark current across the entire detector.

        This assumption is only valid if the scene has a single, centered point 
        source. Any other configuration should trigger a warning.

        """
        # define 2x2 kernel
        kernel = np.asarray([[0.25,0.25],[0.25,0.25]], dtype=np.float32)

        chips = self.current_instrument.get_chip_dimensions()

        brightest_pixel_rate_dict = {}
        for chip in chips:
            temp_rate = rate[chip["bounds"]["y"],chip["bounds"]["x"]]

            # bg_pix_rate contains both bg and dc rates. They SHOULD ultimately both be single-valued
            bg = np.mean(bg_pix_rate[chip["bounds"]["y"],chip["bounds"]["x"]])
            temp_rate += bg

            # convolve - each output value is the sum of 4 adjacent pixels
            brightconv = sg.fftconvolve(temp_rate, kernel, mode='valid') # we only want to pair up pixels in the image; there's no need to extend beyond valid combinations
            # TODO: This is currently dropped on the floor.

            # measure
            bploc = np.unravel_index(np.argmax(temp_rate), temp_rate.shape)
            if self.projection_type == "spec":
                if self.dispersion_axis == "x":
                    wavemax = self.wave_pix[bploc[1]]
                else:
                    wavemax = self.wave_pix[bploc[0]]
            elif self.projection_type == "image":
                wavemax = self.wave[0]
            elif self.projection_type in ("slitless", "multiorder"):
                # The idea of a wavelength-of-brightest-pixel is problematic for slitless and
                # multiorder modes. It can only be reasonably defined for a single centered 
                # single-order source. It fails when a.) the target source is not centered, 
                # b.) the brightest pixel isn't IN the target source (it won't crash if that 
                # brightest pixel is still within the bounds of a centered source, but it'll 
                # be wrong), or c.) when the mode uses CombinedSignal and there's more than 
                # one wavelength grid superimposed on the image.
                # The values produced by this function are not used in JWST or Roman results
                # and as of 2022.0201 we have no HST modes that exercise this particular 
                # clause.
                wavemax = np.nan
            else:
                raise EngineInputError("Unrecognized projection type")
            brightest_pixel_rate_dict[chip["name"]] = (temp_rate[bploc], wavemax)

        return brightest_pixel_rate_dict

    def trigger_saturation_calcs(self):
        # this code must be kept in sync with CombinedSignal.trigger_saturation_calcs()
        # Note that the 2D image due to background alone may have spatial structure due to instrumental effects.
        # Therefore it is calculated here.
        self.bg_pix_rate = self.get_bg_pix_rate()

        # Check to see if the background is saturating (on unbinned data)
        bgsat = self.get_saturation_mask(self.bg_pix_rate)
        if (np.sum(bgsat) > 0) or (np.isnan(np.sum(bgsat))):
            key = "background_saturated"
            self.warnings[key] = warning_messages[key]
        # bin to match the rest of the data
        self.bg_pix_rate = self.binning.sum(self.bg_pix_rate)

        self.ngroup_map, self.sat_type = self.the_detector.get_before_sat(self.rate_plus_bg_unbinned)
        self.ngroup_map = self.binning.min(self.ngroup_map)
        self.fraction_saturation = np.max(self.the_detector.get_saturation_fraction(self.rate_plus_bg_unbinned))

class DetectorBinning(object):
    def __init__(self, rate, bin_x, bin_y, oldgrid, dispersion_axis="x", wave=None, projection_type=""):
        """
        This class exists to model on-detector binning, where individual pixels are exposed
        and then combined with their neighbors within the readout electronics. It can handle
        arbitrary positive integer binning values.

        When instantiated, the class computes the parameters needed to crop and rebin data, 
        as well as the new coordinate IrregularGrid. 
        Binning can then be applied to any 2D detector product that uses the same coordinate grid
        by running one of the methods. Because of the nature of the binning, multiple methods are 
        offered. 
        For instance: We wish to conserve flux within a binned image, so we use the sum() method 
        to sum the pixels.
        However, saturation handling must be different: we can no longer use the fullwell value
        to determine the binning; the fundamental thing that is saturating is the original pixel 
        electronics, and thus a 2x2 binned pixel with 3x the fullwell flux might be constituted 
        of four pixels that each reached 75% fullwell (in which case the flux is not saturated), or 
        it could be a single highly saturated pixel (in which case the flux for that binned pixel is 
        unreliable). Accordingly, we need to compute saturation on the unbinned array,
        and then set each bin to the max() of the saturation of the constituent pixels.

        The basic premise is that binning a 2D image requires three operations:
        1. Crop the array to an even multiple of the binning value
        2. Bin the cropped array
        3. Rewrite the spatial grid so it matches the new cropped, binned dimensions

        Binning a wavelength array requires four operations:
        1. Crop the array to an even multiple of the binning value
        2. Bin the cropped array
        3. Rewrite the spatial grid so it matches the new cropped, binned dimensions
        4. Rewrite the wavelength array using the mean of the values of the binned pixels
        With the added stipulation that for slitless spectroscopy, we are concerned
        with making the slitless spectroscopy trace completely fill the first binned pixel
        (so that we don't have to crop it later)
        
        The geometry of the slitless mode is basically that of the imaging scene with 
        an expanded section in the middle corresponding to the wavelength array: 
        [   *   ] -> [   --------   ]
        We need to make sure that middle section lines up so that the slitless 
        spectrum starts on a pixel: we lose at most one partially-covered pixel at 
        the high-wavelength end that way.

        Methods:
        --------
        sum:
            Rebins and sums the pixels that form each bin. This conserves flux for signal products
        mean:
            Rebins and takes the mean of the pixels that form each bin. This conserves mask 
            shape for extraction masks
        max:
            Rebins and takes the max of the pixels that form each bin (for the saturation mask)
        min:
            Rebins and takes the min of the pixels that form each bin (for the groups-before-sat image)
        dispersion: 
            Rebins the wavelength array, taking the mean of the wavelengths that form each binned pixel.
        """
        # your basic binning needs to be cropped to the nearest multiple of the binning

        self.cropshape = [int(np.floor(rate.shape[0]/bin_y) * bin_y), int(np.floor(rate.shape[1]/bin_x) * bin_x)]
        # and then, of course, binned.
        self.newshape = [int(self.cropshape[0]/bin_y), int(self.cropshape[1]/bin_x)]

        # by default we will assume we want to center the result - trim evenly off the edges
        # y is backwards
        self.end_y = int(np.floor((rate.shape[0]-self.cropshape[0])/2))
        self.start_x = int(np.floor((rate.shape[1]-self.cropshape[1])/2))
        # Now it gets complicated: what kind of rebinning are we doing? image, slit, or slitless?
        if wave is not None:
            if dispersion_axis == 'x':
                bin_wave = bin_x
                self.cropwaveshape = int(np.floor(wave.shape[0]/bin_wave) * bin_wave)
                self.newwaveshape = int(np.floor(wave.shape[0]/bin_wave))
                self.start_wave = self.start_x
                if projection_type=="slitless":
                    # This next equation finds how many pixels there are before the 
                    # wavelength array starts. If that's an even multiple of the binning
                    # we need to do nothing. 
                    self.start_wave = 0
                    offset = int((self.cropshape[1]-len(wave))/2 % bin_wave)
                    if offset % bin_wave > 0:
                        self.start_x += offset
                        # startx was already set to the best position; moving it
                        # most likely means losing another pixel from the far side
                        self.cropshape[1] -= bin_wave
                        self.newshape[1] -= 1 
            else:
                bin_wave = bin_y
                self.cropwaveshape = int(np.floor(wave.shape[0]/bin_wave) * bin_wave)
                self.newwaveshape = int(np.floor(wave.shape[0]/bin_wave))
                self.start_wave = self.end_y
                if projection_type=="slitless":
                    self.start_wave = 0
                    # we want to line up the low-wavelength end, which is at the high-index end
                    offset = bin_wave - int((self.cropshape[0]+len(wave))/2 % bin_wave)
                    if offset % bin_wave > 0:
                        self.end_y += offset
                        self.cropshape[0] -= bin_wave
                        self.newshape[0] -= 1 

        # now, rewrite the Grid instance as an IrregularGrid.
        # Find the actual min x and y position (relating to the edges of the pixels)
        gridy = oldgrid.col[self.end_y: self.end_y+self.cropshape[0]]
        gridx = oldgrid.row[self.start_x: self.start_x+self.cropshape[1]]
        start_y = np.min(gridy) - 0.5 * oldgrid.ysamp
        start_x = np.min(gridx) - 0.5 * oldgrid.xsamp

        # And the new y sampling and number of pixels
        newsampy = bin_y * oldgrid.ysamp
        newsampx = bin_x * oldgrid.xsamp
        # create a new array that starts in the middle of the first pixel and goes to the middle of the last pixel
        gridyvals = np.linspace(start_y + 0.5 * newsampy, start_y + newsampy * (self.newshape[0] - 0.5), self.newshape[0])[::-1]
        gridxvals = np.linspace(start_x + 0.5 * newsampx, start_x + newsampx * (self.newshape[1] - 0.5), self.newshape[1])

        self.grid = coords.IrregularGrid(gridyvals, gridxvals)

    def sum(self, rate):
        # this MUST be sliced in numpy notation ([x,y]) rather than sequentially [x][y]) 
        croprate = rate[self.end_y:self.end_y+self.cropshape[0], self.start_x:self.start_x+self.cropshape[1]]

        # rebin the 2D array
        sh = int(self.newshape[0]), int(self.cropshape[0] // self.newshape[0]), int(self.newshape[1]), int(self.cropshape[1] // self.newshape[1])
        return croprate.reshape(sh).sum(-1).sum(1)

    def mean(self, rate):
        # this MUST be sliced in numpy notation ([x,y]) rather than sequentially [x][y]) 
        croprate = rate[self.end_y:self.end_y+self.cropshape[0], self.start_x:self.start_x+self.cropshape[1]]

        # rebin the 2D array
        sh = int(self.newshape[0]), int(self.cropshape[0] // self.newshape[0]), int(self.newshape[1]), int(self.cropshape[1] // self.newshape[1])
        return croprate.reshape(sh).mean(-1).mean(1)

    def min(self, rate):
        # this MUST be sliced in numpy notation ([x,y]) rather than sequentially [x][y]) 
        croprate = rate[self.end_y:self.end_y+self.cropshape[0], self.start_x:self.start_x+self.cropshape[1]]

        # rebin the 2D array
        sh = int(self.newshape[0]), int(self.cropshape[0] // self.newshape[0]), int(self.newshape[1]), int(self.cropshape[1] // self.newshape[1])
        return croprate.reshape(sh).min(-1).min(1)

    def max(self, rate):
        # this MUST be sliced in numpy notation ([x,y]) rather than sequentially [x][y]) 
        croprate = rate[self.end_y:self.end_y+self.cropshape[0], self.start_x:self.start_x+self.cropshape[1]]

        # rebin the 2D array
        sh = int(self.newshape[0]), int(self.cropshape[0] // self.newshape[0]), int(self.newshape[1]), int(self.cropshape[1] // self.newshape[1])
        return croprate.reshape(sh).max(-1).max(1)

    def dispersion(self, wave):
        if hasattr(self, "start_wave"):
            # We have now set up the cropping required.
            cropwave = wave[self.start_wave:self.start_wave+self.cropwaveshape]
            # rebin the wavelength array
            wavesh = int(self.newwaveshape), int(self.cropwaveshape // self.newwaveshape)
            # we take the mean rather than the sum, because we want the new midpoint
            # wavelength arrays generally aren't linear, so some info is lost here.
            return cropwave.reshape(wavesh).mean(-1)
        else:
            return wave

class CombinedSignal(object):
    """
    This class takes a set of DetectorSignal instances, combines the rates appropriately, and
    provides the other information that DetectorNoise requires.  The primary use for this is in the case
    of SOSS where the detector plane contains signals from effectively three instrument configurations, one
    for each of the visible orders of the gr700xd disperser. These signals need to be combined properly before
    being used to calculate a DetectorNoise and perform a Strategy extraction.

    WARNING: This is currently set up to only support a single aperture slice. Supporting multiple orders with
    multiple slices will require further refactoring...

    Parameters
    ----------
    signal_list: list of DetectorSignal instances
    """

    def __init__(self, signal_list):
        # each signal potentially has a different grid so we want to combine everything onto one grid. instruments
        # that combine multiple signals need to provide detector_pixels so that we know how to shift the signals
        # with respect to each other.
        self.warnings = {}
        self.detector_pixel_list = []
        self.wave_pix_dict = {}
        # some things are common to all signals so get them from the first one
        self.parent_signal = signal_list[0]
        self.signal_list = signal_list
        self.dispersion_axis = self.parent_signal.dispersion_axis
        self.binning = self.parent_signal.binning

        if len(self.parent_signal.rate_list) > 1:
            msg = "Combining multiple instrument signals onto single detector only supports a single slice."
            raise NotImplementedError(value=msg)

        maxx = 0
        maxy = 0
        maxx_unbinned = 0
        maxy_unbinned = 0
        for i, s in enumerate(self.signal_list):
            self.warnings.update(s.warnings)
            pixels = s.detector_pixels
            if pixels is None:
                inst = s.current_instrument.instrument['instrument']
                msg = "No detector pixel configuration set for instrument %s." % inst
                raise DataError(value=msg)
            self.detector_pixel_list.append(pixels)
            ny, nx = s.rate.shape
            maxx = max(maxx, nx)
            maxy = max(maxy, ny)
            ny, nx = s.rate_plus_bg_unbinned.shape
            maxx_unbinned = max(maxx_unbinned, nx)
            maxy_unbinned = max(maxy_unbinned, ny)
            self.wave_pix_dict[s.order] = s.wave_pix

        if self.parent_signal.projection_type == "multiorder":
            self.grid = coords.IrregularGrid(np.arange(maxy), np.arange(maxx) - maxx)
            self.grid_unbinned = coords.IrregularGrid(np.arange(maxy_unbinned), np.arange(maxx_unbinned) - maxx_unbinned)
        else:
            if self.dispersion_axis == 'x':
                # for slitless calculations, the dispersion axis is longer than the spectrum being dispersed
                # because the whole field of view is being dispersed. 'excess' is the size of the FOV
                # and half will be to the left of the blue end of the spectrum and half to the right of the red end.
                # this is used to create the new spatial coordinate transform for the pixel image on the detector.
                self.grid = coords.IrregularGrid(
                    self.parent_signal.grid.col,
                    (np.arange(maxx) - (maxx) / 2.0) * self.parent_signal.grid.xsamp
                )
                self.grid_unbinned = coords.IrregularGrid(
                    self.parent_signal.grid_unbinned.col,
                    (np.arange(maxx_unbinned) - (maxx_unbinned) / 2.0) * self.parent_signal.grid_unbinned.xsamp
                )
            else:
                self.grid = coords.IrregularGrid(
                    (np.arange(maxy) - (maxy) / 2.0) * self.parent_signal.grid.ysamp,
                    self.parent_signal.grid.row
                )
                self.grid_unbinned = coords.IrregularGrid(
                    (np.arange(maxy_unbinned) - (maxy_unbinned) / 2.0) * self.parent_signal.grid_unbinned.ysamp,
                    self.parent_signal.grid_unbinned.row
                )

        self.pixgrid_list = [self.grid]
        self.dist = self.grid.dist()
        self.dist_unbinned = self.grid_unbinned.dist()
        self.det_mask = np.ones_like(self.dist)

        # these parameters are the same for each signal
        self.wave = self.parent_signal.wave
        self.total_flux = self.parent_signal.total_flux
        self.fp_rate = self.parent_signal.fp_rate
        self.bg_fp_rate = self.parent_signal.bg_fp_rate
        self.background = self.parent_signal.background
        self.flux_cube_list = self.parent_signal.flux_cube_list
        self.flux_plus_bg_list = self.parent_signal.flux_plus_bg_list

        self.aperture_list = self.parent_signal.aperture_list
        self.projection_type = self.parent_signal.projection_type
        self.current_instrument = self.parent_signal.current_instrument
        self.spatial_grid = self.parent_signal.grid
        self.spatial_grid_unbinned = self.parent_signal.grid_unbinned
        self.the_detector = self.parent_signal.the_detector
        self.calculation_config = self.parent_signal.calculation_config
        if self.parent_signal.the_detector.rn_correlation:
            self.read_noise_correlation_matrix = self.parent_signal.read_noise_correlation_matrix

        self.unbinned_array = np.zeros_like(self.dist_unbinned)

        self.rate_list = [{
            'fp_pix': np.zeros_like(self.dist),
            'fp_pix_no_ipc': np.zeros_like(self.dist),
            'fp_pix_no_qyield': np.zeros_like(self.dist),
            'fp_pix_variance': np.zeros_like(self.dist),
            'fp_pix_unbinned': np.zeros_like(self.dist_unbinned),
            'fp_pix_no_ipc_unbinned': np.zeros_like(self.dist_unbinned),
            'fp_pix_no_qyield_unbinned': np.zeros_like(self.dist_unbinned),
            'fp_pix_variance_unbinned': np.zeros_like(self.dist_unbinned)
        }]
        self.rate_plus_bg_list = [{
            'fp_pix': np.zeros_like(self.dist),
            'fp_pix_no_ipc': np.zeros_like(self.dist),
            'fp_pix_no_qyield': np.zeros_like(self.dist),
            'fp_pix_variance': np.zeros_like(self.dist),
            'fp_pix_unbinned': np.zeros_like(self.dist_unbinned),
            'fp_pix_no_ipc_unbinned': np.zeros_like(self.dist_unbinned),
            'fp_pix_no_qyield_unbinned': np.zeros_like(self.dist_unbinned),
            'fp_pix_variance_unbinned': np.zeros_like(self.dist_unbinned)
        }]
        saturation = np.zeros_like(self.dist_unbinned)

        for i, s in enumerate(self.signal_list):
            # each fp_pix and fp_pix_no_ipc product already has background and dark current added.
            # to avoid double-counting the dark current, we need to subtract it from each piece
            # and add it back in once at the end.
            s.rate_plus_bg_list[0]['fp_pix_unbinned'] -= s.the_detector.dark_current
            s.rate_plus_bg_list[0]['fp_pix_no_ipc_unbinned'] -= s.the_detector.dark_current
            s.rate_plus_bg_list[0]['fp_pix_variance_unbinned'] -= s.the_detector.dark_current

            ny, nx = s.rate_list[0]['fp_pix_unbinned'].shape
            # set up position within the combined rate image to add each rate
            if self.dispersion_axis == 'x':
                ly = 0
                uy = int(self.unbinned_array.shape[0] - s.rate_list[0]['fp_pix_unbinned'].shape[0])
                # shifts between the rates are referenced to the first pixel location
                lx = 0
                ux = int(self.unbinned_array.shape[1] - s.rate_list[0]['fp_pix_unbinned'].shape[1])
            else:
                # shifts between the rates are referenced to the first pixel location
                ly = 0
                uy = int(self.unbinned_array.shape[0] - s.rate_list[0]['fp_pix_unbinned'].shape[0])
                lx = 0
                ux = int(self.unbinned_array.shape[1] - s.rate_list[0]['fp_pix_unbinned'].shape[1])

            # this will cause pixels that don't overlap to be filled with the background. this is correct for imaging
            # and slitless, but may need revisiting for other cases.
            new_r = np.pad(s.rate_list[0]['fp_pix_unbinned'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_bg = np.pad(s.rate_plus_bg_list[0]['fp_pix_unbinned'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_noipc = np.pad(s.rate_list[0]['fp_pix_no_ipc_unbinned'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_bg_noipc = np.pad(s.rate_plus_bg_list[0]['fp_pix_no_ipc_unbinned'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_noqyield = np.pad(s.rate_list[0]['fp_pix_no_qyield_unbinned'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_bg_noqyield = np.pad(s.rate_plus_bg_list[0]['fp_pix_no_qyield_unbinned'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_var = np.pad(s.rate_list[0]['fp_pix_variance_unbinned'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_bg_var = np.pad(s.rate_plus_bg_list[0]['fp_pix_variance_unbinned'], ([ly, uy], [lx, ux]), mode='edge')
            sat_mask = np.zeros_like(s.rate_plus_bg_list[0]['fp_pix_no_ipc_unbinned'])
            saturation_mask = self.the_detector._get_saturation_mask(s.rate_plus_bg_list[0]['fp_pix_no_ipc_unbinned'], sat_mask)
            new_sat = np.pad(saturation_mask, ([ly, uy], [lx, ux]), mode='constant')

            saturation = np.maximum(saturation, new_sat)

            self.rate_list[0]['fp_pix_unbinned'] += new_r
            self.rate_list[0]['fp_pix_no_ipc_unbinned'] += new_r_noipc
            self.rate_list[0]['fp_pix_no_qyield_unbinned'] += new_r_noqyield
            self.rate_list[0]['fp_pix_variance_unbinned'] += new_r_var
            self.rate_plus_bg_list[0]['fp_pix_unbinned'] += new_r_bg
            self.rate_plus_bg_list[0]['fp_pix_no_ipc_unbinned'] += new_r_bg_noipc
            self.rate_plus_bg_list[0]['fp_pix_no_qyield_unbinned'] += new_r_bg_noqyield
            self.rate_plus_bg_list[0]['fp_pix_variance_unbinned'] += new_r_bg_var

        # now we add the dark current back in
        # this assumes a scalar dark level; there is no allowance for a non-uniform dark to be added
        # back here.
        self.rate_plus_bg_list[0]['fp_pix_unbinned'] += self.parent_signal.the_detector.dark_current
        self.rate_plus_bg_list[0]['fp_pix_no_ipc_unbinned'] += self.parent_signal.the_detector.dark_current
        self.rate_plus_bg_list[0]['fp_pix_variance_unbinned'] += self.parent_signal.the_detector.dark_current

        self.rate_list[0]['fp_pix'] = self.binning.sum(self.rate_list[0]['fp_pix_unbinned'])
        self.rate_list[0]['fp_pix_no_ipc'] = self.binning.sum(self.rate_list[0]['fp_pix_no_ipc_unbinned'])
        self.rate_list[0]['fp_pix_no_qyield'] = self.binning.sum(self.rate_list[0]['fp_pix_no_qyield_unbinned'])
        self.rate_list[0]['fp_pix_variance'] = self.binning.sum(self.rate_list[0]['fp_pix_variance_unbinned'])
        self.rate_plus_bg_list[0]['fp_pix'] = self.binning.sum(self.rate_plus_bg_list[0]['fp_pix_unbinned'])
        self.rate_plus_bg_list[0]['fp_pix_no_ipc'] = self.binning.sum(self.rate_plus_bg_list[0]['fp_pix_no_ipc_unbinned'])
        self.rate_plus_bg_list[0]['fp_pix_no_qyield'] = self.binning.sum(self.rate_plus_bg_list[0]['fp_pix_no_qyield_unbinned'])
        self.rate_plus_bg_list[0]['fp_pix_variance'] = self.binning.sum(self.rate_plus_bg_list[0]['fp_pix_variance_unbinned'])

        self.saturation_list = [self.binning.max(saturation)]
        # SOSS only has one aperture so the on_detector rates are just the fp_pix rates.
        self.rate = self.rate_list[0]['fp_pix']
        self.rate_no_qyield = self.rate_list[0]['fp_pix_no_qyield']
        self.rate_plus_bg = self.rate_plus_bg_list[0]['fp_pix']
        self.rate_unbinned = self.rate_list[0]['fp_pix_unbinned']
        self.rate_plus_bg_unbinned = self.rate_plus_bg_list[0]['fp_pix_unbinned']

        self.trigger_saturation_calcs()
        self.brightest_pixel = np.max(self.rate_plus_bg)
        self.groups_list = [self.ngroup_map]

        # HST Health and Safety
        self.brightest_pixel_rate = self.parent_signal.calc_brightest_pixel(self.rate_no_qyield,self.bg_pix_rate)
        self.detector_total_rate = self.parent_signal.calc_total_detector(self.rate_no_qyield,self.bg_pix_rate)

    def get_saturation_mask(self, rate=None):
        """
        Compute a numpy array indicating pixels with full saturation (2), partial saturation (1) and no saturation (0).
        This version just wraps whats implemented within DetectorSignal.

        Parameters
        ----------
        rate: None or 2D np.ndarray
            Detector plane rate image used to build saturation map from

        Returns
        -------
        mask: 2D np.ndarray
            Saturation mask image
        """
        if rate is None:
            rate = self.rate_plus_bg_unbinned

        mask = self.parent_signal.get_saturation_mask(rate=rate)
        return mask

    def spectral_detector_transform(self):
        """
        Create engine API format dict section containing properties of wavelength coordinates
        at the detector plane.

        Returns
        -------
        t: dict (engine API compliant keys)
        """
        return self.parent_signal.spectral_detector_transform()

    def spectral_model_transform(self):
        """
        Create engine API format dict section containing properties of the wavelength coordinates
        used in the construction of a ConvolvedSceneCube.

        Returns
        -------
        t: dict (engine API compliant keys)
        """
        return self.parent_signal.spectral_model_transform()

    def cube_wcs_info(self):
        """
        Create WCS headers and FITS binary table that describe the cube's coordinate system.

        Returns
        -------
        tbhdu: astropy.io.fits.BinTableHDU instance
            The wavelength sampling is irregularly spaced so we define a binary FITS
            table that contains the array of wavelengths, self.wave.
        header: dict
            Contains the WCS keys that define the coordinate transformation for all axes
        """
        return self.parent_signal.cube_wcs_info()

    def trigger_saturation_calcs(self):
        # this code must be kept in sync with DetectorSignal.trigger_saturation_calcs()
        self.bg_pix_rate = self.rate_plus_bg_unbinned - self.rate_unbinned
        # Check to see if the background is saturating (on unbinned data)
        bgsat = self.get_saturation_mask(self.bg_pix_rate)
        if (np.sum(bgsat) > 0) or (np.isnan(np.sum(bgsat))):
            key = "background_saturated"
            self.warnings[key] = warning_messages[key]
        # bin to match the rest of the data
        self.bg_pix_rate = self.binning.sum(self.bg_pix_rate)

        # create ngroup map from the newly combined signal
        self.ngroup_map, self.sat_type = self.signal_list[0].the_detector.get_before_sat(self.rate_plus_bg_unbinned)
        self.ngroup_map = self.binning.min(self.ngroup_map)
        self.fraction_saturation = np.max(
            self.signal_list[0].the_detector.get_saturation_fraction(self.rate_plus_bg_unbinned))


    def export_to_fits(self, fitsfile='ModelDetectorCube'):
        """
        Write convolved scene cube to a FITS file
        """
        header = self.grid.wcs_info()
        slice_index = 0
        for s,signal in enumerate(self.signal_list):
            for slice_index,flux_plus_bg in enumerate(signal.flux_plus_bg_list):
                fitsfile_slice = '{}_{}_{}.fits'.format(fitsfile, str(slice_index).strip(), s)
                t1 = fits.PrimaryHDU(np.moveaxis(flux_plus_bg, [0, 1, 2], [1, 2, 0]))
                c1 = fits.Column(name='wavelength', array=signal.wave, format='D')
                t2 = fits.BinTableHDU.from_columns([c1])
                t1.header.update(header)
                tbhdu = fits.HDUList([t1,t2])
                tbhdu.writeto(fitsfile_slice, overwrite=True)
                tbhdu.close()
