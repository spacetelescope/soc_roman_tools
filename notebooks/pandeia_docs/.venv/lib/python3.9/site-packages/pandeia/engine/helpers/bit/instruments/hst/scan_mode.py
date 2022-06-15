""" This module contains functions specific to scan mode.

"""

from __future__ import division

from ...engine import initialize
from ...engine.area_over_which import ImageRectangularArea
from ...engine.string_constants import SNR_REGION, SPECTROSCOPIC, SCSPECTROSCOPIC, IMAGING, SCIMAGING


scan_science_modes = {
    IMAGING: SCIMAGING,
    SPECTROSCOPIC: SCSPECTROSCOPIC
}


def make_rectangular_img_snr_region(args, instrument):
    """ Makes a specialized imaging extraction region
        used in scan mode.

    Parameters
    ----------
    args: tuple or list
      with the region dimesnions and units
    instrument: an instance of Imager
      the instrument instance

    Returns
    -------
    instance of ImageRectangularArea

    """
    width = args[1][0]
    height = args[1][1]
    units = args[2]

    arcsec_per_pixel_height = instrument.detector.pixel_angheight
    arcsec_per_pixel_width = instrument.detector.pixel_angwidth

    rgn = ImageRectangularArea(SNR_REGION,
                               width,
                               height,
                               units,
                               arcsec_per_pixel_width,
                               arcsec_per_pixel_height)
    return rgn


def get_scan_configuration_data(ei):
    """ Loads instrument and detector data from config files.

        This function calls the base function and then loads
        scan mode parameters into the custom configuration dict.

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
    iconfig, dconfig, custom_config = initialize.get_configuration_data(ei)

    # scan mode parameters
    custom_config['scan_mode'] = ei['science_mode'] in (SCIMAGING, SCSPECTROSCOPIC)
    if custom_config['scan_mode']:
        custom_config['scan_length'] = ei['sclength']

    return iconfig, dconfig, custom_config


def modify_scan_configuration(instrument, custom_config):
    """ The scan_mode flag must be propagated into the instrument
        instance. Science mode must be translated properly, since
        the instruments by default only recognize 'imaging' and
        'spectroscopic' as valid modes. As a consequence, the
        limits config file must be revisited again.

    Parameters
    ----------
    instrument: Imager or Spectrograph
      Instance being configured
    custom_config: dictionary
      custom configuration info

    """
    if custom_config['scan_mode']:

        instrument.scan_mode = True
        instrument.scan_length = custom_config['scan_length']
        instrument.science_mode = scan_science_modes[instrument.science_mode]

        # In scan mode, the limits config file must be revisited because
        # the science mode was just redefined by statement above. The
        # scan mode-specific limits may have been missed by the default
        # read that happened in the generic instrument._configure method.
        instrument.set_instrument_detector_limits()
