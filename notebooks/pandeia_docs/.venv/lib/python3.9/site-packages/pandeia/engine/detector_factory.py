# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

from .custom_exceptions import EngineInputError
from .utils import recursive_subclasses, merge_data

from .detector import Detector, Detector_MultiAccum, Detector_SingleAccum
from .detector import H2RG, H4RG, SiAs
from .detector import CCD, MAMA, H1R, XDL

def DetectorFactory(webapp=False, config={}, **kwargs):
    """
    Function to take a detector type and configuration data and build/return an appropriately
    configured instance of the desired Detector class.

    Parameters
    ----------
    webapp: bool
        Switch to toggle strict API checking
    config: dict
        Configuration data segments needed to set up the exposure calculation
    **kwargs: keyword/value pairs
        Additional configuration data
    """

    types = recursive_subclasses(Detector)
    detectors = [t.__name__.lower() for t in types]
    type_map = dict(list(zip(detectors, types)))

    exposure_string = config["detector_config"]["det_type"].lower()

    cls = type_map[exposure_string](webapp=webapp, config=config, **kwargs)
    return cls
