"""
Function signatures:
--------------------

``compatibility(engine_input, result)``
    No return value. 'result' is the engine_input after compatibility changes are applied.

``get_configuration_data(ei)``
    Returns instrument configuration dictionaries given the engine input dictionary.

``get_detector_tables(instrument_name, detector_name, science_mode, spec_element)``
    Returns a dictionary with the enclosed energy tables associated with the disperser.

``modify_areas_of_interest(instrument)``
    No return value; modifies instrument.areas_of_interest

``make_snr_region(instrument, *args, **kwargs)`` #different for imager or spectrograph
    Returns an AreaOverWhich of the correct type; modifies instrument.areas_of_interest

``make_area_pixel(instrument, *args, **kwargs)`` #different for imager or spectrograph?
    Returns an AreaOverWhich for 1-pix extraction; modifies instrument.areas_of_interest

``modify_configuration(instrument, custom_config)``
    No return value. Modifies the Instrument instance in place.

``modify_target_or_sky(instrument, observed_target, observed_sky, obswave=None)``
    Returns new instances of ObservedItem with the modified target and
    sky, or references to the input, if no modification took place.

``stray_light(instrument, target, sky)``
    Returns a ObservedItem with the stray light component, or None.

``make_readpattern(instrument, engine_input)``
    Returns a ReadPattern instance configured for the instrument and user inputs.

``modify_noise(noise_collection)``
    No return value. Modifies the NoiseCollection instance in place.

``modify_total(instrument, total_collection)``
    No return value. Modifies the total DetectedCollection instance in place.

``modify_one_pixel(instrument, extracted_collection)``
    Returns a modified copy of the input one-pixel ExtractedCollection.

``detect_dark(instrument, dark)``
    Returns a DetectedItem or DetectedCollection with the value of the dark current per pixel.

``prepare_vectors_for_plotting(observed_target, background)``
    Modifies target and/or background vector objects according to
    the expectations of what the plots and table should do.
"""
from __future__ import division

from .telescopes import TELESCOPES
#from ..main.tools import find_instrument
from pandeia.engine.helpers.bit.pyetc_util import dynamic_import, stable_dict, find_instrument

MODULE_PATH_NAME = 'pandeia.engine.helpers.bit.instruments.%s.%s.instr_functions'

# Custom configuration functions handle instrument-specific configuration data
#
ALLOWED_CUSTOM_CONFIG = ['compatibility',
                         'get_configuration_data',
                         'modify_configuration'
                        ]

# Custom calculation functions implement instrument-specific calculations.
#
# Although the detector tables are read as configuration information,
# they actually represent PSF calculations that are expressed via table
# interpolation. Hence their inclusion in the subset of calculation functions.
ALLOWED_CUSTOM_CALC = ['get_detector_tables',
                       'make_snr_region',
                       'modify_areas_of_interest',
                       'modify_target_or_sky',
                       'stray_light',
                       'make_area_whole_detector',
                       'make_area_pixel',
                       '_init_custom_config',
                       'compute_buffer_time',
                       'make_readpattern',
                       'modify_total',
                       'modify_noise',
                       'modify_one_pixel',
                       'detect_dark',
                       'prepare_vectors_for_plotting']

def import_custom_overrides(instrument_name):
    """ Imports and makes available to the caller a dictionary of
        instrument-specific functions.

    Parameters
    ----------

    instrument_name: str
       name of instrument

    Returns
    -------

    result: dict
        dict with instrument-specific functions.

    """
    result = stable_dict()

#    try:
#        instr_module = find_instrument(TELESCOPES[instrument_name], instrument_name + ".instr_functions")

        # workaround for the above
    module_name = MODULE_PATH_NAME % (TELESCOPES[instrument_name], instrument_name)
    instr_module = dynamic_import(module_name)
#    except ImportError, e:
#        # No custom imports
#        return result

    for key in ALLOWED_CUSTOM_CONFIG:
        try:
            result[key] = instr_module.__getattribute__(key)
        except AttributeError as e:
            continue

    for key in ALLOWED_CUSTOM_CALC:
        try:
            result[key] = instr_module.__getattribute__(key)
        except AttributeError as e:
            continue

    return result

