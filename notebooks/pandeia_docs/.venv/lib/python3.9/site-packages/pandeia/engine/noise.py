# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import numpy as np
from warnings import warn

from .custom_exceptions import DataError

class DetectorNoise(object):
    """
    This class provides the functionality to model the noise in a JWST observation

    Parameters
    ----------
    obs_signal: DetectorSignal or CombinedSignal instance
        The calculated signal at the detector plane
    observation: observation.Observation instance
        Configuration information for the Observation that generated ObsSignal
    """

    def __init__(self, obs_signal, observation):
        self.warnings = {}
        self.observation = observation

        self.grid = obs_signal.grid
        self.dist = obs_signal.dist
        self.det_mask = obs_signal.det_mask
        self.calculation_config = obs_signal.calculation_config
        self.binning = obs_signal.binning

        self.current_instrument = observation.instrument
        self.the_detector = self.current_instrument.the_detector

        # Are we dealing with an image or spectroscopic mode?
        self.projection_type = self.current_instrument.projection_type

        # get saturation so we can calculate unsaturated groups
        # if (a if condition else b)
        # if the calculation_config value is True or False (edited by user), use True or False. If it's none, use the instrument-team-defined default from detector_config
        if (self.calculation_config.effects['saturation'] if self.calculation_config.effects['saturation'] is not None else self.the_detector.saturation):
            self.the_detector.set_saturation(True)
        else:
            self.the_detector.set_saturation(False)

        # if (a if condition else b)
        # if the calculation_config value is True or False (edited by user), use True or False. If it's none, use the instrument-team-defined default from detector_config
        if (self.calculation_config.noise['crs'] if self.calculation_config.noise['crs'] is not None else self.the_detector.crs):
            # if we are including CRs, we need a CR rate configured for the telescope. this rate is in events/s/cm^2/sr.
            if not hasattr(self.current_instrument.telescope, 'cr_rate'):
                msg = "No cosmic ray rate defined for %s." % self.current_instrument.telescope.tel_name.upper()
                raise DataError(value=msg)
            self.cr_rate = self.current_instrument.telescope.cr_rate

            # we need physical pixel size to calculate incident area in sq cm
            if not hasattr(self.the_detector, 'pix_size_x') or not hasattr(self.the_detector, 'pix_size_y'):
                msg = "The physical pixel size of the detector in microns must be specified to calculate the CR rate."
                raise DataError(value=msg)
            self.pix_size_x = self.the_detector.pix_size_x
            self.pix_size_y = self.the_detector.pix_size_y

            # need to scale the rate by number of pixels affected by each event.
            if not hasattr(self.current_instrument.telescope, 'cr_npixels'):
                msg = "No defined value for number of pixels affected per CR event for %s." % \
                      self.current_instrument.telescope.tel_name.upper()
                raise DataError(value=msg)
            self.cr_npixels = self.current_instrument.telescope.cr_npixels

            cr_rate = self.cr_rate * self.cr_npixels

            # cr_rate is in units of events/s/cm^2/sr.  convert pix_size in um to cm and multiply by 4pi steradians
            # to get cr_rate per pixel.
            self.pix_cr_rate = (self.pix_size_x * 1.0e-4 * self.pix_size_y * 1.0e-4) * 4.0 * np.pi * cr_rate
        else:
            self.pix_cr_rate = 0.0

        self.var_pix_list, self.stdev_pix_list, self.rn_var_pix_list, self.ff_pix_list = self.basic_source_noise(obs_signal)
        self.var_pix, self.stdev_pix, self.var_rn_pix, self.ff_pix = self.on_detector()

    def on_detector(self):
        """
        Calculates the detector plane noise products

        Returns
        -------
        products: tuple
            detector_var - numpy.ndarray
                Full variance including all noise sources
            detector_stdev - numpy.ndarray
                Standard deviation (i.e. sqrt(detector_var))
            detector_rn_var - numpy.ndarray
                Variance stricly due to detector readnoise
            detector_ff_var - numpy ndarray
                Variance strictly due to flat field noise
        """
        aperture_sh = self.var_pix_list[0].shape
        n_apertures = len(self.var_pix_list)
        detector_shape = (aperture_sh[0] * n_apertures, aperture_sh[1])
        detector_var = np.zeros(detector_shape)
        detector_stdev = np.zeros(detector_shape)
        detector_rn_var = np.zeros(detector_shape)
        detector_ff = np.zeros(detector_shape)

        i = 0
        for var_pix, stdev_pix, rn_var_pix, ff_pix in zip(self.var_pix_list, self.stdev_pix_list, self.rn_var_pix_list, self.ff_pix_list):
            detector_var[i * aperture_sh[0]:(i + 1) * aperture_sh[0]:] = var_pix
            detector_stdev[i * aperture_sh[0]:(i + 1) * aperture_sh[0]:] = stdev_pix
            detector_rn_var[i * aperture_sh[0]:(i + 1) * aperture_sh[0]:] = rn_var_pix
            detector_ff[i * aperture_sh[0]:(i + 1) * aperture_sh[0]:] = ff_pix

            i += 1

        output = self.det_mask * detector_var, self.det_mask * detector_stdev, self.det_mask * detector_rn_var, self.det_mask * detector_ff
        return output

    def basic_source_noise(self, obs_signal):
        """
        Calculate the noise using the full pixelated flux cube.

        Inputs
        ------
        obs_signal: ObsSignal class instance
            Class containing the detector flux plane or plane cube

        Returns
        -------
        var_pix_list: List of ndarrays
            The detector variance plane or plane cube
        stdev_pix_list: List of ndarrays
            The detector standard deviation plane or plane cube (literally the square root of var_pix_list).
        rn_var_pix_list: List of ndarrays
            Variance strictly due to detector readnoise
        ffnoise_pix_list: List of ndarrays
            Variance strictly due to flat field noise
        """
        ff_electrons = self.the_detector.ff_electrons

        var_pix_list = []
        stdev_pix_list = []
        rn_var_pix_list = []
        ffnoise_pix_list = []

        for rate_plus_bg in obs_signal.rate_plus_bg_list:
            slope_var_pix, slope_rn_var_pix = self.get_slope_variance(rate_plus_bg)
            rate_per_pix = rate_plus_bg['fp_pix']
            """
            The flat field error is a division by ~1 (the flat field is normalized), with a variance of 1/ff_electrons.
            Note that the value of the flat field response is constant for multiple ramps and multiple integrations, so
            nramps > 1 does not decrease the residual flat field noise. Due to that, a user will either have to improve the
            flat field or dither with > 1 pixel offsets. The most apparent effect for everyday ETC use is that this sets an
            upper limit on the achievable signal-to-noise ratio.

            The pixel variance upon division with a normalized flat field constructed with ff_electrons (it is assumed that
            the flat field is ideal):

            s^2(R/FF) = s^2(R) + R^2/FF_electrons
            """
            var_pix = slope_var_pix / self.the_detector.exposure_spec.nramps
            rn_var_pix = slope_rn_var_pix / self.the_detector.exposure_spec.nramps

            # Add the flat field residual noise if requested. If there are multiple exposures, decrease the noise
            # to account for the effect of co-adding exposures.
            # if (a if condition else b)
            # if the calculation_config value is True or False (edited by user), use True or False. If it's none, use the instrument-team-defined default from detector_config
            if (self.calculation_config.noise['ffnoise'] if self.calculation_config.noise['ffnoise'] is not None else self.the_detector.ffnoise):
                ffnoise_pix = (rate_per_pix ** 2 / ff_electrons) / self.the_detector.exposure_spec.nexp
                var_pix += ffnoise_pix
            else:
                ffnoise_pix = np.zeros_like(rate_per_pix)

            stdev_pix = np.sqrt(var_pix)

            var_pix_list.append(var_pix)
            stdev_pix_list.append(stdev_pix)
            rn_var_pix_list.append(rn_var_pix)
            ffnoise_pix_list.append(ffnoise_pix)

        products = var_pix_list, stdev_pix_list, rn_var_pix_list, ffnoise_pix_list
        return products


    def get_slope_variance(self, rate):
        """
        Compute the slope variance, based on detector properties and
        exposure specifications.

        Inputs
        ------
        rate_per_pix: ndarray
           The measured slope of the ramp (ie, the rate) per pixel
           This can be a one-dimensional wavelength spectrum (for a 1D ETC)
        or a three-dimensional cube ([wave,x,y] for a 3D ETC).

        Returns
        -------
        slope_var: ndarray
           The associated variance of the measured slope
        slope_rn_var: ndarray
           The associated variance of the readnoise only
        """

        exp_pars = self.the_detector.exposure_spec


        unsat_ngroups = self.the_detector.get_unsaturated(rate['fp_pix_no_ipc'], full_saturation=self.the_detector.mingroups)
        # scale the unsaturated ngroups by the CR loss. While the fundamental thing being hit is individual unbinned pixels,
        # the math we've already done computes the effect per area (of unbinned pixel) and is thus still accurate
        #
        # The fact that unsat_ngroups is not an integer is intentional, and is the way the ETC statistically
        # accounts for CR losses.
        # if (a if condition else b)
        # if the calculation_config value is True or False (edited by user), use True or False. If it's none, use the instrument-team-defined default from detector_config
        if (self.calculation_config.noise['crs'] if self.calculation_config.noise['crs'] is not None else self.the_detector.crs):
            unsat_ngroups = np.vectorize(self.the_detector.calc_cr_loss)(unsat_ngroups, self.pix_cr_rate)


        slope_var, slope_rn_var = self.the_detector.get_slope_variance(rate, unsat_ngroups)

        return slope_var, slope_rn_var

    def get_readnoise_slope_variance(self, rate):
        """
        Compute the variance of the read noise only, excluding all other noise contributions.

        Inputs
        ------
        none

        Returns
        -------
        var_rn: ndarray
            The variance of the readnoise given the current detector slope parameters. Note that
            if the readnoise is switched off, the associated variance will of course be zero.
        """

        unsat_ngroups = self.the_detector.get_unsaturated(rate['fp_pix_no_ipc'],
                                                        full_saturation=self.mingroups)
        # scale the unsaturated ngroups by the CR loss.
        # if (a if condition else b)
        # if the calculation_config value is True or False (edited by user), use True or False. If it's none, use the instrument-team-defined default from detector_config
        if (self.calculation_config.noise['crs'] if self.calculation_config.noise['crs'] is not None else self.the_detector.crs):
            unsat_ngroups = np.vectorize(self.the_detector.calc_cr_loss)(unsat_ngroups, self.pix_cr_rate)

        var_rn = self.the_detector.rn_variance(self.the_detector.rn, unsat_ngroups=unsat_ngroups)

        return var_rn
