from __future__ import division


""" ACS-specific functions

"""



from ....engine import initialize


def get_configuration_data(ei):
    """ Load instrument and detector data from config files.

        This overrides the base function in order to handle
        ACS specifics.

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

    # spectral resolution must be returned in the detector config dict.
    dconfig['spectral_resolution'] = iconfig[dconfig['detector_name']]['spectral_resolution']

    return iconfig, dconfig, custom_config

