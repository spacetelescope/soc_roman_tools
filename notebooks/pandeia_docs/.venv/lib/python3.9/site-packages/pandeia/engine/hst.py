# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

from .telescope import Telescope
from .instrument import Instrument
from .custom_exceptions import DataError
from .constants import DEFAULT_DISPERSION
import astropy.units as u
import stsynphot as st
import numpy as np
from numpy.polynomial import polynomial


class HST(Telescope):

    def get_ote_eff(self, wave):
        """
        Temporary function; allow stsynphot to handle all of the graph table functionality in one go.
        """
        return np.ones_like(wave)



class HSTInstrument(Instrument):

    """
    Generic HST Instrument class
    """
    def __init__(self, mode=None, config={}, **kwargs):
        telescope = HST()
        self.instrument_pars = {}
        # TODO HST has in general single accum detectors for which ngroup and nint do not make sense.
        self.instrument_pars['detector'] = ["nsplit"]
        self.instrument_pars['instrument'] = ["aperture", "disperser", "detector", "filter", "instrument", "mode"]
        self.api_parameters = list(self.instrument_pars.keys())

        # these are required for calculation, but ok to live with config file defaults
        self.api_ignore = ['dynamic_scene', 'max_scene_size', 'scene_size']
        Instrument.__init__(self, telescope=telescope, mode=mode, config=config, **kwargs)
        # avoid doing the complex obsmode calculations multiple times
        self.obsmode = self.get_obsmode()
        self.rpower, self.dlds, self.dispwave = self.get_dispersion_products()

    def get_obsmode(self):
        if "imaging" in self.instrument['mode']:
            obsmode = self._obsmode_imaging()
        elif "spectroscopy" in self.instrument['mode']:
            obsmode = self._obsmode_spec()
        elif "echelle" in self.instrument['mode']:
            obsmode = self._obsmode_spec()
        elif "target_acq" in self.instrument['mode']:
            obsmode = self._obsmode_ta()
        if "mjd" in self.instrument:
            mjd = self.instrument['mjd']
        else:
            mjd = self.telescope.mjd
        obsmode += ",mjd#{}".format(mjd)

        return obsmode

    def _obsmode_imaging(self):
        """
        Temporary function for basic imaging obsmodes, to be replaced as we work modes.
        These are likely not correct and not complete, but they should be functional.
        """
        return "{},{},{}".format(self.instrument['instrument'], self.instrument['detector'], self.instrument['filter'])

    def _obsmode_spec(self):
        """
        Temporary function for basic spectroscopic obsmodes, to be replaced as we work modes.
        These are likely not correct and not complete, but they should be functional.
        """
        if "cenwave" in self.instrument:
            return "{},{},{},{}".format(self.instrument['instrument'], self.instrument['detector'], self.instrument['disperser'], self.instrument['cenwave'])
        else:
            return "{},{},{}".format(self.instrument['instrument'], self.instrument['detector'], self.instrument['disperser'])

    def get_dispersion_products(self):
        """
        For HST we don't have separate files, we have the TRDS graph table
        and the default wavelength range being output
        """

        band = st.band(self.obsmode)
        dispwave = band.binset.to_value(u.micron)

        # code from pandeia_data/devtools/create_MIRI_LRS_dispersion.py
        dlds_int = np.gradient(dispwave)
        dlds_coeffs = polynomial.polyfit(dispwave,dlds_int,7)
        dlds = polynomial.polyval(dispwave, dlds_coeffs)

        r = dispwave/dlds
        return r, dlds, dispwave

    def get_dispersion(self, wave):
        """
        Return dispersion

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate/trim dispersion data onto

        Returns
        -------
        dispersion: numpy.ndarray or float
            Dispersion as a function of wave
        """
        disperser = self._get_disperser()
        if disperser is not None:
            dlds = self.dlds
            dlds = self._interp_refdata(wave, dlds, "dlds", default=DEFAULT_DISPERSION)
        else:
            dlds = 0.0
        return dlds

    def get_resolving_power(self, wave):
        """
        Return the resolving power

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate resolving power onto

        Returns
        -------
        rpower: numpy.ndarray or None
            Resolving power as a function of wave
        """
        disperser = self._get_disperser()
        if disperser is not None:
            rpower = self.rpower
            rpower = self._interp_refdata(wave, rpower, "R", default=wave/(DEFAULT_DISPERSION*2))
        else:
            rpower = None
        return rpower

    def get_wave_pix(self):
        """
        Get the wavelength vector to convert pixel position to wavelength
        For HST, stsynphot and the graph table have already done this for us.

        Returns
        -------
        wavepix: numpy.ndarray
            Wavelength vector mapping pixels to wavelengths
        """
        return self.dispwave

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
            pixels = np.arange(len(wave))
        else:
            pixels = None
        return pixels

    def get_wave_blaze(self):
        """
        Return the wavelength vector used in the grating efficiency (blaze) file

        Returns
        -------
        wave_blaze: numpy.ndarray
            Wavelength vector from the grating efficiency file
        """
        return self.dispwave


    def _interp_refdata(self, wave, data, name, default=None):
        """
        Read in requested reference file and return requested data interpolated onto wave

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength array onto which the reference data is to be interpolated
        data: numpy.ndarray
            Key to get filename from self.paths
        name: string
            Name of the reference data; only used for the error message
        default: None or float
            Default value to return in case of missing file or column

        Returns
        -------
        ref: numpy.ndarray or float
            If file exists, return reference_data(wave). Else return default.
        """
        ref = None
        if data is not None:
            ref = np.interp(wave, self.dispwave, data)
        if ref is None and default is None:
            msg = "No reference data found and no default value supplied for %s." % (name)
            raise DataError(value=msg)
        elif ref is None and default is not None:
            ref = np.ones_like(wave)*default
        return ref

    def _get_cdbs(self, wave, key):
        """
        HST data is available in the $PYSYN_CBDS trees, so it should be loaded differently.

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength array onto which the throughput curve is to be interpolated
        key: str
            Key to get filename from the $PYSYN_CDBS tree

        Returns
        -------
        eff: numpy.ndarray or float
            If ref file exists, return efficiency(wave), else return 1.0
        """
        if 'None' not in key:
            bp = st.spectrum.band(key)
            eff = bp(wave*u.micron).value #don't want to pass the Quantity
        else:
            # If it's explicitly None, there is no element; pass through 1.0
            eff = np.ones_like(wave)*1.0
        return eff

    def get_filter_eff(self, wave):
        """
        Temporary function; allow stsynphot to handle all of the graph table functionality in one go.
        """
        return np.ones_like(wave)

    def get_disperser_eff(self, wave):
        """
        Temporary function; allow stsynphot to handle all of the graph table functionality in one go.
        """
        return np.ones_like(wave)

    def get_detector_qe(self, wave):
        """
        Temporary function; allow stsynphot to handle all of the graph table functionality in one go.
        """
        return np.ones_like(wave)

    def get_internal_eff(self,wave):
        """
        This is the only temporary function that will actually do anything. It relies on
        stsynphot for graph table functionality
        
        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        eff: numpy.ndarray or float
            Disperser efficiency as a function of wave
        """
        obsmode = self.get_obsmode()
        bp = st.spectrum.band(obsmode)
        eff = bp(wave*u.micron).value #don't want to pass the Quantity

        return eff

    def get_quantum_yield(self,wave):
        """
        Compute detector quantum yield if the critical wavelength is defined. If not, 
        return unity (trivial yield).

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        quantum_yield: numpy.ndarray or float
            Quantum efficiency as a function of wave
        """
        if hasattr(self.the_detector, "critical_wavelength"):
            q_yield = self.the_detector.critical_wavelength / wave
        else:
            q_yield = 1
        fano_factor = (3.0 * q_yield - q_yield**2. - 2.) / q_yield

        return q_yield, fano_factor

class COS(HSTInstrument):

    """
    Special methods unique to HST COS
    """

    def _loadpsfs(self):
        """
        The COS PSFs differ greatly by focus position, and we thus need to 
        know lifetime position (currently only LP4), disperser, and cenwave.
        """
        if self.instrument['mode'] == "fuv_spectroscopy":
            psf_key = "{}{}{}{}".format(self.instrument['aperture'], "lp4", self.instrument['disperser'], self.instrument['cenwave'])
        elif self.instrument['mode'] == "nuv_imaging":
            psf_key = self.instrument['aperture']+self.instrument['slit']
        self.psf_library = self._load_psf_library(psf_key)

    def _obsmode_imaging(self):
        items = ['cos',self.instrument['aperture'],self.instrument['detector'], self.instrument['slit']]
        obsmode = ",".join(items)
        return obsmode


class WFC3(HSTInstrument):

    """
    Special methods unique to HST WFC3
    """
    pass

class ACS(HSTInstrument):

    """
    Special methods unique to HST ACS
    """
    pass

class STIS(HSTInstrument):

    """
    Special methods unique to HST STIS
    """
    def __init__(self, mode=None, config={}, **kwargs):
        """
        STIS has additional configuration parameters, like fuv_glow_region and dark_level.
        Do the standard HSTInstrument setup, then read in the extra parameters.
        """
        super().__init__(mode, config, **kwargs)

        if 'fuvmama' in self.instrument['detector']:
            fuv_glow_key = 'dark_glow_region_{}'.format(config['detector']['fuv_glow_region'])
            fuv_glow_rate = self.detector_config['fuvmama'][fuv_glow_key]
            self.the_detector.dark_current = self.the_detector.dark_current + fuv_glow_rate

    def _obsmode_imaging(self):
        items = ['stis',self.instrument['detector'], 'mirror']
        if self.instrument['filter'] is not None:
            items.append(self.instrument['filter'])
        obsmode = ",".join(items)
        return obsmode

    def _obsmode_spec(self):
        """
        Much of the PyETC version of this code has to deal with the fact that the web interface
        does not contain the prefix character (i or c). In pandeia, the names of the cenwaves include
        the prefix.
        """
        disperser = self.instrument['disperser']
        detector = self.instrument['detector']
        obsmode = ",".join(['stis', detector, disperser])

        cenwave = self.instrument.get('cenwave',None)
        if cenwave != None:
            obsmode += ",{}".format(cenwave)
        else: #use default central wavelength if there is one
            try:
                cenwave = self.config_constraints['disperser'][disperser]['cenwaves'][0]
                obsmode += ",{}".format(cenwave)
            except (KeyError, IndexError):
                pass
                #otherwise leave it off
        
        slit = self.instrument.get('slit', None)
        if slit is not None:
            obsmode += ",{}".format(slit)

        return obsmode

    def apply_scattering(self, extracted_list):
        """
        Compute and apply the echelle scattered light parameter
        Will not affect any configuration that does not have
        defined echelle scattering 

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto
        rate: extracted target rate

        """
        if hasattr(self.the_detector, "echelle_scattering"):
            wave_scatter, mult_scatter = self.the_detector.echelle_scattering
            for idx,extracted in enumerate(extracted_list):

                scattering_factor = np.interp(extracted["wavelength"], wave_scatter, mult_scatter)

                scatter_rate = (1. / scattering_factor - 1.) * extracted["extracted_noise"]

                # at this point, the noise is sqrt(variance), but we need to add this 
                # as a new poisson noise term.
                extracted_list[idx]["extracted_noise"] = np.sqrt(extracted["extracted_noise"]**2 + scatter_rate)

        return extracted_list

    def get_global_scattering(self):
        """
        Return the global scattering rate, scattered portion
        as defined in the data
        """
        cenwave = self.instrument.get('cenwave', None)

        # No cenwave or cenwave = all will give no detector scattering, but PyETC 
        # doesn't make that option possible.
        scatter_rate = 0
        if hasattr(self.the_detector, "echelle_global"):
            if cenwave in self.the_detector.echelle_global:
                scatter_rate = self.the_detector.echelle_global[cenwave] - 1

        return scatter_rate
