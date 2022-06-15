# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import numpy as np
import numpy.ma as ma
from .custom_exceptions import EngineInputError, DataError

class ExposureSpec:
    """
    Parent class for all exposure times, MultiAccum and SingleAccum
    """
class ExposureSpec_MultiAccum(ExposureSpec):
    """
    Parent class for MultiAccum
    """

    def __init__(self, config={}, webapp=False, **kwargs):
        """
        Create a generic Exposure Specification.

        Inputs
        ------
        config: dict
            dictionary of detector configuration setups

        webapp: bool
            Switch to toggle strict API checking
            
        **kwargs: keyword/value pairs
            Additional configuration data
        """
        self.webapp = webapp

        # Required parameters
        self.readout_pattern = config["input_detector"]["readout_pattern"]
        self.subarray = config["input_detector"]["subarray"]
        self.ngroup = config["input_detector"]["ngroup"]
        self.nint = config["input_detector"]["nint"]
        self.nexp = config["input_detector"]["nexp"]
        self.tframe = config["subarray"]["default"][self.subarray]["tframe"]
        self.tfffr = config["subarray"]["default"][self.subarray]["tfffr"]
        if self.readout_pattern in config["subarray"] and self.subarray in config["subarray"][self.readout_pattern]:
            self.tframe = config["subarray"][self.readout_pattern][self.subarray]["tframe"]
            self.tfffr = config["subarray"][self.readout_pattern][self.subarray]["tfffr"]

        # Optional parameters
        # These are defined by the Instrument's reference data as Instrument properties.
        self.nframe = config["readout_pattern"][self.readout_pattern]["nframe"]
        self.ndrop2 = config["readout_pattern"][self.readout_pattern]["ndrop2"]

        self.nprerej = config["detector_config"]["nprerej"]
        self.npostrej = config["detector_config"]["npostrej"]
        # Going from general to specific: start with the default reset values, then check the defaults, then check
        # the actual readout pattern for this subarray. Use the most appropriate reset values for this setup, otherwise
        # the defaults
        self.nreset1 = 1
        self.nreset2 = 1
        if "nreset1" in config["subarray"]["default"][self.subarray]:
            self.nreset1 = config["subarray"]["default"][self.subarray]["nreset1"]
            self.nreset2 = config["subarray"]["default"][self.subarray]["nreset2"]
        if self.readout_pattern in config["subarray"] and self.subarray in config["subarray"][self.readout_pattern]:
            if "nreset1" in config["subarray"][self.readout_pattern][self.subarray]:
                self.nreset1 = config["subarray"][self.readout_pattern][self.subarray]["nreset1"]
                self.nreset2 = config["subarray"][self.readout_pattern][self.subarray]["nreset2"]
        if "nreset1" in config["readout_pattern"][self.readout_pattern]:
            self.nreset1 = config["readout_pattern"][self.readout_pattern]["nreset1"]
            self.nreset2 = config["readout_pattern"][self.readout_pattern]["nreset2"]


        # These are never specified in our data, currently; they were always the default
        # value from the ExposureSpec __init__ signature.
        if "ndrop1" in config["readout_pattern"][self.readout_pattern]:
            self.ndrop1 = config["readout_pattern"][self.readout_pattern]["ndrop1"]
        else:
            self.ndrop1 = 0
        if "ndrop3" in config["readout_pattern"][self.readout_pattern]:
            self.ndrop3 = config["readout_pattern"][self.readout_pattern]["ndrop3"]
        else:
            self.ndrop3 = 0
        self.frame0 = False
        if "frame0" in config["subarray"]["default"][self.subarray]:
            self.frame0 = config["subarray"]["default"][self.subarray]["frame0"]
        if self.readout_pattern in config["subarray"] and self.subarray in config["subarray"][self.readout_pattern]:
            if "frame0" in config["subarray"][self.readout_pattern][self.subarray]:
                self.nreset1 = config["subarray"][self.readout_pattern][self.subarray]["frame0"]

        # If these are trivial, we don't have to define them.
        if "nsample" in config["readout_pattern"][self.readout_pattern]:
            self.nsample = config["readout_pattern"][self.readout_pattern]["nsample"]
            self.nsample_skip = config["readout_pattern"][self.readout_pattern]["nsample_skip"]
        else:
            self.nsample = 1
            self.nsample_skip = 0

        # Target acqs only use a subset of the groups
        if "ngroup_extract" in config["input_detector"]:
            self.ngroup_extract = config["input_detector"]["ngroup_extract"]

        self.get_times()

        # Derived quantities
        self.nramps = self.nint * self.nexp

    def get_times(self):
        """
        The time formulae are defined in Holler et al. 2021, JWST-STScI-006013-A. Note that we have to subtract the
        groups that are rejected (by the pipeline) from the measurement time (nprerej+npostrej). The saturation time
        conservatively considers the ramp saturated even if saturation occurs in a rejected frame.

        Also note that these equations are generic, suitable for both H2RG and SiAs detectors.

        The equations in this method are duplicated in the front-end
        (workbook.js, update_detector_time_labels function). Please ensure
        that changes to these equations are reflected there.
        """

        self.tgroup = self.tframe * (self.nframe + self.ndrop2)

        # MIRI measurement time for ngroups < 5 is now handled with dropframes
        # This reduces to Equation 3 for H2RG detectors and Equation 5 for SiAs detectors.
        self.measurement_time = self.nint * self.tframe * (self.ngroup - 1 - self.nprerej - self.npostrej) * (self.nframe + self.ndrop2)
        # Equation 4
        if self.frame0:
            self.measurement_time += 0.5 * self.nint * self.tframe * (self.nframe - 1)

        # Equation 1, which naturally simplifies to Equation 2 for SiAs detectors.
        self.exposure_time = (self.tfffr * self.nint) + self.tframe * (self.nreset1 + (self.nint - 1) * self.nreset2 +
                                self.nint * (self.ndrop1 + (self.ngroup - 1) * (self.nframe + self.ndrop2) +
                                self.nframe + self.ndrop3))

        # Equation 6, which reduces to Equation 7 for SiAs detectors.
        self.saturation_time = self.tframe * (self.ndrop1 + (self.ngroup - 1) * (self.nframe + self.ndrop2) + self.nframe)

        self.duty_cycle = self.measurement_time / self.exposure_time
        self.total_exposure_time = self.nexp * self.exposure_time

        self.total_integrations = self.nexp * self.nint

class ExposureSpec_H2RG(ExposureSpec_MultiAccum):

    pass


class ExposureSpec_H4RG(ExposureSpec_MultiAccum):

    pass

class ExposureSpec_SiAs(ExposureSpec_MultiAccum):

    def get_times(self):
        super().get_times()

        # This is where we adjust values so we can still use the same
        # MULTIACCUM formula as for the NIR detectors. We need the effective
        # "average time per sample" for MIRI.
        self.tsample = self.tframe / (self.nsample + self.nsample_skip)
        # 'nsample_total' for MIRI is the total number of non-skipped samples X number of averaged frames.
        # Note that in practice, it currently never happens that both the number of samples and number of
        # averaged frames are >1 (i.e., no SLOWGRPAVG exists). However, this will deal with that situation,
        # should it occur.
        self.nsample_total = self.nframe * self.nsample

class ExposureSpec_SingleAccum(ExposureSpec):
    """
    Parent class for SingleAccum
    """

    def __init__(self, config={}, webapp=False, **kwargs):
        """
        Create a single accum Exposure Specification.

        Inputs
        ------
        config: dict
            dictionary of detector configuration setups

        webapp: bool
            Switch to toggle strict API checking

        **kwargs: keyword/value pairs
            Additional configuration data

        """
        self.webapp = webapp

        self.time = config["input_detector"]["time"]

        # Required parameters
        #self.readout_pattern = config["input_detector"]["readout_pattern"]
        #self.subarray = config["input_detector"]["subarray"]
        if "nexp" in config["input_detector"]:
            raise DataError("SingleAccum calculations cannot use nexp")
        self.nsplit = config["input_detector"]["nsplit"]

        self.get_times()

        # Derived quantities needed for the generic noise equation
        self.nramps = 1

    def get_times(self):
        self.measurement_time = self.time
        self.exposure_time = self.time / self.nsplit
        self.saturation_time = self.time / self.nsplit
        self.total_exposure_time = self.time

        self.duty_cycle = 1.
        self.total_integrations = 1


class ExposureSpec_CCD(ExposureSpec_SingleAccum):

    pass

class ExposureSpec_H1R(ExposureSpec_CCD):

    pass

class ExposureSpec_MAMA(ExposureSpec_SingleAccum):

    pass

class ExposureSpec_XDL(ExposureSpec_MAMA):

    pass
