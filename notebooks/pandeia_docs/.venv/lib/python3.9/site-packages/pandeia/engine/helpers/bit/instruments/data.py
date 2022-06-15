from __future__ import division

import os
import importlib # Requires Python 2.7 or higher

import pandeia.engine.helpers.bit.pyetc_util as util

class Data(object):
    """
    Interface class to instrument/detector configuration files

    Data class provides the interface between the instrument/detector
    configuration files and the rest of pyetc.

    Parameters
    ----------
    location: str
        location of data files
    instrument: str
        name of instrument
    telescope: str
        name of telescope

    """
    def __init__(self, telescope, instrument, location=__file__):

        self.location = location
        self.location = "/".join(os.path.abspath(location).split('/')[0:-1])

        self.instrument = instrument
        self.telescope = telescope
        self.dconfig = util.stable_dict()
        self.iconfig = util.stable_dict()
        self.custom_config = util.stable_dict()
        self.aperture_fraction = util.stable_dict()
        self.sharpness = util.stable_dict()
        self.extraction_defaults = util.stable_dict()

        self.populate_iconfig()
        self.populate_dconfig()
        self.populate_custom_config()
        self.populate_aperture_fraction()
        self.populate_sharpness()
        self.populate_extraction_defaults()

        # Extract the relevant limits
        self._inst_limits =  importlib.import_module('pandeia.engine.helpers.bit.instruments.' +
                                               telescope + '.' + instrument +
                                               '.' + instrument + '_limits')

    def populate_iconfig(self):
        """
        Populate the iconfig dictionary for use in engine.initialize

        """
        #Read in the instrument config file
        iconfig_name = self.location +'/' + self.instrument + '_config.dat'
        self.iconfig = util.read_dict(iconfig_name)

    def populate_dconfig(self):
        """
        Populate the dconfig dictionary for use in engine.initialize

        """
        dconfig_name = self.location + '/' + 'detector_config.dat'
        self.dconfig = util.read_dict(dconfig_name)

    def populate_custom_config(self):
        """
        Populate the custom_config dictionary for potential use in engine.initialize

        """

        # A custom configuration file is not required. Check for existence
        # and return an empty dictionary if not found.
        try:
            custom_config_name = self.location + '/' + 'custom_config.dat'
            self.custom_config = util.read_dict(custom_config_name)
        except IOError as e:
            self.custom_config = util.stable_dict()

    def populate_aperture_fraction(self):
        """
        Populate the instrument specific aperture fraction table

        """
        #Read in the aperture fraction data file
        try:
            aperture_fraction_name = self.location +'/' + self.instrument + '_aperture_fraction.dat'
            self.aperture_fraction = util.read_dict(aperture_fraction_name)
        except IOError as e:
            self.aperture_fraction = util.stable_dict()

    def populate_sharpness(self):
        """
        Populate the instrument specific sharpness table

        """
        #Read in the sharpness data file
        try:
            sharpness_name = self.location +'/' + self.instrument + '_sharpness.dat'
            self.sharpness = util.read_dict(sharpness_name)
        except IOError as e:
            self.sharpness = util.stable_dict()

    def populate_extraction_defaults(self):
        """
        Populate the extraction default dictionary if available

        """
        #Read in the extraction defaults data file
        try:
            extraction_defaults_name = self.location +'/' + self.instrument + '_extraction_defaults.dat'
            self.extraction_defaults = util.read_dict(extraction_defaults_name)
        except IOError as e:
            self.extraction_defaults = util.stable_dict()

    def get_dconfig(self, detector_name):
        """
        Return a copy of the dconfig dictionary

        """
        return self.dconfig[detector_name].copy()

    def get_iconfig(self):
        """
        Return a copy of the iconfig dictionary

        """
        return self.iconfig.copy()

    def get_custom_config(self):
        """
        Return a copy of the custom_config dictionary

        """
        return self.custom_config.copy()

    def get_aperture_fraction(self):
        """
        Return a copy of the aperture_fraction dictionary

        """
        return self.aperture_fraction.copy()

    def get_sharpness(self):
        """
        Return a copy of the sharpness dictionary

        """
        return self.sharpness.copy()

    def get_extraction_defaults(self):
        """
        Return a copy of the extraction_defaults dictionary

        """
        return self.extraction_defaults.copy()

    def get_relevant_limits(self, mode, detector, filter=''):
        """
        Use the engine.limit_parser interface to
        extract the instrument specific limit information

        """
        from ..engine import limit_parser
        limitdict = limit_parser.extract_limits(mode,
                                                detector,
                                                filter,
                                                self._inst_limits.limit_spec,
                                                self._inst_limits.allsub)

        return limitdict
