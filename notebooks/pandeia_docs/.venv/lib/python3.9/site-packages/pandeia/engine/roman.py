# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

from .telescope import Telescope
from .instrument import Instrument

class Roman(Telescope):

    """
    Currently a dummy class for directory/file discovery, but could eventually contain Roman-specific methods
    """
    pass


class RomanInstrument(Instrument):

    """
    Generic Roman Instrument class
    """

    def __init__(self, mode=None, config={}, **kwargs):
        telescope = Roman()
        # these are the required sections and need to be passed via API in webapp mode
        self.instrument_pars = {}
        self.instrument_pars['detector'] = ["nexp", "ngroup", "nint", "readout_pattern", "subarray"]
        self.instrument_pars['instrument'] = ["aperture", "disperser", "filter", "instrument", "mode"]
        self.api_parameters = list(self.instrument_pars.keys())

        # these are required for calculation, but ok to live with config file defaults
        self.api_ignore = ['dynamic_scene', 'max_scene_size', 'scene_size']

        Instrument.__init__(self, telescope=telescope, mode=mode, config=config, **kwargs)

    def read_detector(self):
        """Read in the detector keyword from the aperture parameter in the config json file.
           Put the detector keyword in self.instrument['detector']."""
        self.instrument['detector'] = self.aperture_config[self.get_aperture()]['detector']


class WFI(RomanInstrument):

    """
    Currently, the Roman WFI requires only one method beyond those provided by the generic Instrument class
    """
    def _loadpsfs(self):
        """
        Short-wavelength filters need PSFs with the skinny pupil mask.
        Long-wavelength filters (just f184 and f213 for now) need the wide pupil mask.
        The grism and prism both have their own pupil masks.
        Because the ranges overlap, we need to switch between PSF libraries.
        The mask is specifically baked into the filters, so selecting by filter seems appropriate
        """
        if self.instrument['mode'] == 'imaging':
            psf_key = self.instrument['filter']
        elif self.instrument['mode'] == 'spectroscopy':
            psf_key = self.instrument['disperser']
        else:
            message = "Invalid mode specification: {}".format(str(self.instrument['mode']))
            raise EngineInputError(value=message)

        self.psf_library = self._load_psf_library(psf_key)


class IFU(RomanInstrument):

    """
    Currently the Roman IFU (deprecated) requires no extra methods beyond those provided by the generic Instrument class
    """
    pass
