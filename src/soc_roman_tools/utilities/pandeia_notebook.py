"""
Some utilites to work with Pandeia.engine and Roman.

    Written by Russell Ryan

    NOTES.... things changed for v1.7.1
              1) the aperture is now iehter imaging or spectroscopy
              2) subarray is now are imaging, defocus_small, defocus_large,
                 for imaging or spectroscopy for spectroscopy

    calc['configuration']['instrument']['aperture'] = 'imaging'
    calc['configuration']['instrument']['disperser'] = None
    calc['configuration']['detector']['readout_pattern'] = readout_pattern
    calc['configuration']['detector']['subarray'] = 'imaging'
"""

from scipy.optimize import minimize_scalar
import pkg_resources
import os
import warnings

# Pandeia stuff
from pandeia.engine.perform_calculation import perform_calculation
from pandeia.engine.calc_utils import build_default_calc

# Require this version of Pandeia
PANDEIA_REQ = '1.7'
VALID_FILTERS = ('f062', 'f087', 'f106', 'f129', 'f146', 'f158', 'f184', 'f213')

DEFAULTS = {'geometry': 'point', 'sed_type': 'flat', 'unit': 'fnu',
            'norm_type': 'at_lambda', 'norm_fluxunit': 'abmag', 'norm_wave': 2.0,
            'norm_waveunit': 'um', 'readout_pattern': 'medium8', 'ngroup': 10,
            'nexp': 1, 'nint': 1,
            'background_level': 'low', 'background': None,  # (0.9,0.),
            'aperture_size': 0.2, 'sky_annulus': (0.6, 0.8),
            'background_subtraction': True}


def set_parameter(kwargs, param):
    """
    Purpose
    -------
    Simple function to use pass value or DEFAULT value.
    """

    if param in kwargs:
        return kwargs[param]
    elif param in DEFAULTS:
        return DEFAULTS[param]
    else:
        raise ValueError(f"keyword {param} not found")


def pandeia_version():
    """
    Purpose
    -------
    Get the version of the Pandeia code and reference data.
    """

    # Pandeia has this functionality, but it only prints to STDOUT.  So,
    # I'm following the same steps to get the code and data versions
    # pandeia.engine.pandeia_version()

    versions = {'code': pkg_resources.get_distribution('pandeia.engine').version}

    # Get the data version
    with open("{}/VERSION_PSF".format(os.environ['pandeia_refdata'])) as fp:
        # last char is a newline
        versions['data'] = fp.readline().strip()

    return versions


def make_scene(mag=25., **kwargs):
    """
    Purpose
    -------
    Make a default scene dictionary.

    Inputs
    ------
    mag (float):
        Magnitude of the flux normalization.

    Returns
    -------
    scene (dictionary):
        A dictionary containing the scene parameters for pandeia.
    """

    scene = {'shape': {'geometry': set_parameter(kwargs, 'geometry')},
             'spectrum': {'spectrum_parameters': ['sed', 'normalization']}}
    scene['spectrum']['sed'] = {'sed_type': set_parameter(kwargs, 'sed_type'),
                                'unit': set_parameter(kwargs, 'unit')}
    scene['spectrum']['normalization'] = {'type': set_parameter(kwargs, 'norm_type')}
    scene['spectrum']['normalization']['norm_flux'] = mag
    scene['spectrum']['normalization']['norm_fluxunit'] = \
        set_parameter(kwargs, 'norm_fluxunit')
    scene['spectrum']['normalization']['norm_wave'] = \
        set_parameter(kwargs, 'norm_wave')
    scene['spectrum']['normalization']['norm_waveunit'] = \
        set_parameter(kwargs, 'norm_waveunit')

    return scene


def make_calculation(filt, **kwargs):
    """
    Purpose
    -------
    Make a default calculation dictionary.

    Inputs
    ------
    filt (string):
        Roman optical element name.

    Returns
    -------
    calc (dictionary):
        A dictionary containing the calculation parameters for pandeia.
    """
    filt = filt.lower()
    assert (filt in VALID_FILTERS), 'Invalid filter entry'

    # It's Roman, so start with that.
    calc = build_default_calc('roman', 'wfi', 'imaging')

    # make a default scene into this calculation
    calc['scene'][0] = make_scene(**kwargs)

    # some a bevy of things, some of which shouldn't be altered
    calc['configuration']['instrument']['filter'] = filt
    calc['configuration']['instrument']['aperture'] = 'imaging'
    calc['configuration']['instrument']['disperser'] = None
    calc['configuration']['detector']['readout_pattern'] = \
        set_parameter(kwargs, 'readout_pattern')
    calc['configuration']['detector']['subarray'] = 'imaging'
    calc['configuration']['detector']['ngroup'] = set_parameter(kwargs, 'ngroup')
    calc['configuration']['detector']['nint'] = set_parameter(kwargs, 'nint')
    calc['configuration']['detector']['nexp'] = set_parameter(kwargs, 'nexp')

    # set background levels
    calc['background_level'] = set_parameter(kwargs, 'background_level')
    calc['strategy']['target_xy'] = [0.0, 0.0]
    calc['strategy']['aperture_size'] = set_parameter(kwargs, 'aperture_size')
    calc['strategy']['sky_annulus'] = set_parameter(kwargs, 'sky_annulus')
    calc['strategy']['units'] = 'arcsec'  # THIS MUST ALWAYS BE ARCSEC
    calc['strategy']['background'] = set_parameter(kwargs, 'background')
    calc['strategy']['background_subtraction'] = \
        set_parameter(kwargs, 'background_subtraction')

    return calc


def _metric_(sntarg, snmeas):
    """
    Purpose
    -------
    Helper function. This is to hone in on a "target" S/N.

    Inputs
    ------
    sntarg (float):
        Target signal-to-noise.

    snmeas (float):
        Measured signal-to-noise.

    Returns
    -------
    square_diff (float):
        Square of the difference of the target and measured signal-to-noise
        ratios.
    """

    square_diff = (sntarg - snmeas) ** 2

    return square_diff


def compute_sn(filt, mag, nexp, **kwargs):
    """
    Purpose
    -------
    Method to compute the signal-to-noise ratio given a magnitude and number
    of exposures.

    Inputs
    ------
    filt (string):
        Roman optical element name.

    mag (float):
        Flux normalization magnitude.

    nexp (integer):
        Number of exposures.

    Returns
    -------
    sn (float):
        Signal-to-noise value from the ETC calculation.

    etc (~pandeia.engine.report.Report):
        Result of the pandeia ETC calculation.
    """

    calc = make_calculation(filt, **kwargs, mag=mag, nexp=nexp)

    etc = perform_calculation(calc)
    sn = etc['scalar']['sn']

    return sn, etc


def _mag2sn_(mag, calc, sntarget):
    """
    Purpose
    -------
    Helper function to optimize the signal-to-noise given a magnitude.

    Inputs
    ------
    mag (float):
        Flux normalization magnitude.

    calc (dictionary):
        A dictionary containing the pandeia calculation parameters.

    sntarget (float):
        Target signal-to-noise ratio.

    Returns
    -------
    metric (float):
        The squared difference between the measured and target signal-to-noise
        ratios.
    """

    calc['scene'][0]['spectrum']['normalization']['norm_flux'] = mag
    etc = perform_calculation(calc)['scalar']
    metric = _metric_(sntarget, etc['sn'])

    return metric


def compute_mag(filt, sn, nexp, bracket=(18, 30), xtol=1e-4, **kwargs):
    """
    Purpose
    -------
    Compute the magnitude from the signal-to-noise ratio and number of exposures.

    Inputs
    ------
    filt (string):
        Roman optical element name.

    sn (float):
        Signal-to-noise ratio.

    nexp (integer):
        Number of exposures.

    bracket (tuple of floats; optional; default=(18, 30)):
        A tuple containing the lower and upper bounds for the magnitude calculation
        result.

    xtol (float; optional; default=1e-4):
        Tolerance for the magnitude calculation result.

    Returns
    -------
    mag (float):
        Magnitude result from the ETC calculation.

    etc (~pandeia.engine.report.Report):
        Result of the Pandeia ETC calculation.
    """

    calc = make_calculation(filt, nexp=nexp, **kwargs)

    res = minimize_scalar(_mag2sn_, bracket=bracket, bounds=bracket, args=(calc, sn),
                          method='brent', options={'xtol': xtol})
    mag = res['x']
    calc['scene'][0]['spectrum']['normalization']['norm_flux'] = mag
    etc = perform_calculation(calc)

    return mag, etc


def _nexp2sn_(nexp, calc, sntarget):
    """
    Purpose
    -------
    Helper function to minimize the S/N as function of the number
    of exposures.

    Inputs
    ------
    nexp (integer):
        Number of exposures.

    calc (dictionary):
        A dictionary containing the Pandeia calculation parameters.

    sntarget (float):
        Target signal-to-noise ratio.

    Returns
    -------
    metric (float):
        The squared difference between the measured and target signal-to-noise
        ratios.
    """

    calc['configuration']['detector']['nexp'] = int(nexp)
    etc = perform_calculation(calc)['scalar']
    metric = _metric_(sntarget, etc['sn'])

    return metric


def compute_nexp(filt, sn, mag, bracket=(1, 1000), xtol=0.1, **kwargs):
    """
    Purpose
    -------
    Method to compute the number of exposures from a requested signal-to-noise
    ratio and source magnitude.

    Inputs
    ------
    filt (string):
        Roman optical element name.

    sn (float):
        Signal-to-noise ratio.

    mag (flaot):
        Flux normalization magnitude.

    bracket (tuple of floats; optional; default=(1, 1000)):
        A tuple containing the lower and upper bounds for the number of exposures
        calculation result.

    xtol (float; optional; default=0.1):
        Tolerance for the number of exposures calculation result.

    Returns
    -------
    nexp (integer):
        Number of exposures.

    etc (~pandeia.engine.report.Report):
        Result of the Pandeia ETC calculation.
    """

    calc = make_calculation(filt, mag=mag, **kwargs)

    # just check the minimum (if we get a S/N higher than requested in
    # a single exposure, then there's no sense in trying to minimize
    # things, as 1 exposure is the shortest quantum of exposure time
    # (since we don't fluctate Ngroup)
    calc['configuration']['detector']['nexp'] = 1
    etc = perform_calculation(calc)
    if etc['scalar']['sn'] > sn:

        warnings.warn('The S/N for a single exposure is *LARGER* than requested. '
                      'Therefore returning nexp=1 and the S/N is *more* than '
                      'requested.')

        nexp = 1
    else:

        res = minimize_scalar(_nexp2sn_, bracket=bracket, bounds=bracket,
                              args=(calc, sn), method='bounded',
                              options={'xatol': xtol})
        nexp = int(res['x'])
        calc['configuration']['detector']['nexp'] = nexp
        etc = perform_calculation(calc)

        # this generally returns a S/N less than the required amount.
        # let's ensure that we get *AT LEAST* the required S/N for 2 reasons:
        # 1) better to err on the side of caution
        # 2) make code consistent with the above if-clause
        if etc['scalar']['sn'] < sn:
            nexp += 1

    return nexp, etc
