from __future__ import division


""" WFC3IR-specific functions

"""



from .. import scan_mode

from ....engine.imager import Imager
from ....engine.spectrograph import Spectrograph
from ....engine.string_constants import RECTANGLE, SNR_REGION, PIX_REGION
from ....engine.area_over_which import SpecPixAOW


def make_snr_region(instrument, *args, **kwargs):
    """
    WFC3IR in scan mode uses specialized SNR and
    1-pix regions.

    """
    rgn = None

    if isinstance(instrument, Imager) and instrument.scan_mode:

        rgn = scan_mode.make_rectangular_img_snr_region(args, instrument)

    elif isinstance(instrument, Spectrograph) and instrument.scan_mode:

        height = args[0]
        units = args[1]
        obswave = args[2]

        if obswave is None:
            raise ValueError("Must specify obswave.")

        rgn = SpecPixAOW(SNR_REGION, height, units, obswave,
                         instrument._pixels_per_resel,
                         instrument.detector.pixel_angheight,
                         instrument.detector.pixel_angwidth)
        pix_rgn = SpecPixAOW(PIX_REGION, height, units, obswave,
                             1,
                             instrument.detector.pixel_angheight,
                             instrument.detector.pixel_angwidth)

        instrument.areas_of_interest[PIX_REGION] = pix_rgn

    return rgn


def get_configuration_data(ei):
    """ Load instrument and detector data from config files.

        This overrides the base function in order to handle
        WFC3IR specifics.

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
