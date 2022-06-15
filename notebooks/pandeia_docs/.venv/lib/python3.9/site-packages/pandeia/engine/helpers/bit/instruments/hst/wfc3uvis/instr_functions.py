from __future__ import division


""" WFC3UVIS-specific functions

"""



from .. import scan_mode

from ....engine.imager import Imager
from ....engine.string_constants import RECTANGLE


def make_snr_region(instrument, *args, **kwargs):
    """
    WFC3UVIS imaging in scan mode uses a specialized SNR
    rectangular region.

    """
    if instrument.scan_mode:
        return scan_mode.make_rectangular_img_snr_region(args, instrument)
    else:
        return None


def get_configuration_data(ei):
    """ Load instrument and detector data from config files.

        This overrides the base function in order to handle
        WFC3UVIS specifics.

    Parameters
    ----------
    ei: dict
      engine input

    Returns
    -------
    iconfig: dict
      dictionary of instrument data
    dconfig: dict
      dictionary of detector data
    custom_config: dict
      dictionary with custom data

    """
    # get baseline configuration data.
    return scan_mode.get_scan_configuration_data(ei)


def modify_configuration(instrument, custom_config):
    """ The scan_mode flag must be propagated into the instrument
        instance, and also from the raw table data into instances
        of EE table.

    Parameters
    ----------
    instrument: Imager or Spectrograph
      Instance being configured
    custom_config: dictionary
      custom configuration info

    """
    scan_mode.modify_scan_configuration(instrument, custom_config)
