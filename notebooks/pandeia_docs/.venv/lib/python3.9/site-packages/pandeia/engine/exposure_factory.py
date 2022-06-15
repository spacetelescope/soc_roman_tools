# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

from .custom_exceptions import EngineInputError
from .utils import recursive_subclasses, merge_data

# need to import the Instrument subclasses from whereever they're defined
from .exposure import ExposureSpec, ExposureSpec_MultiAccum, ExposureSpec_SingleAccum
from .exposure import ExposureSpec_H2RG, ExposureSpec_H4RG, ExposureSpec_SiAs
from .exposure import ExposureSpec_CCD, ExposureSpec_H1R, ExposureSpec_MAMA, ExposureSpec_XDL

def ExposureFactory(det_type, webapp=False, config={}, **kwargs):
    """
    Function to take a detector type and configuration data and build/return an appropriately
    configured instance of the desired Exposure class.

    Parameters
    ----------
    det_type: str
        The type of detector we want to make exposures with
    webapp: bool
        Switch to toggle strict API checking
    config: dict
        Configuration data segments needed to set up the exposure calculation
    **kwargs: keyword/value pairs
        Additional configuration data
    """
    method = None
    all_config = merge_data(config, dict(**kwargs))
    types = recursive_subclasses(ExposureSpec)
    detectors = [t.__name__.lower() for t in types]
    type_map = dict(list(zip(detectors, types)))

    exposure_string = "exposurespec_{}".format(det_type.lower())

    cls = type_map[exposure_string](webapp=webapp, config=config, **kwargs)
    return cls
