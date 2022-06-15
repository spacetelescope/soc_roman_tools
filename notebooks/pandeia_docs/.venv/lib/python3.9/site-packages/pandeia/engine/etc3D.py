# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import copy
import numpy as np
from warnings import warn

from . import observation
from .report import ReportFactory
from .scene import Scene
from .calc_utils import build_empty_scene
from .custom_exceptions import EngineInputError
from .instrument_factory import InstrumentFactory
from .strategy import StrategyFactory
from .pandeia_warnings import etc3d_warning_messages as warning_messages
from . import debug_utils
from .signal import DetectorSignal, CombinedSignal, CalculationConfig
from .noise import DetectorNoise


# pyfftw.interfaces.cache.enable()

def setup(input, webapp=False):
    """
    This generic portion of the code unpacks the input dictionary and uses the 
    pieces to configure the class instances that will perform the calculation.
    """

    warnings = {}
    try:
        scene_configuration = input['scene']
        background = input['background']
        background_level = input.get('background_level', 'medium')
        instrument_configuration = input['configuration']
        strategy_configuration = input['strategy']
        if input.get('debugarrays'):
            debug_utils.init(input.get('debugarrays'))
    except KeyError as e:
        message = "Missing information required for the calculation: %s" % str(e)
        raise EngineInputError(value=message)

    # get the calculation configuration from the input or use the defaults
    if 'calculation' in input:
        calc_config = CalculationConfig(config=input['calculation'])
    else:
        calc_config = CalculationConfig()

    # #### BEGIN calculation #### #
    """
    This section currently implements the Pandeia engine API.  As the engine's object model
    is refactored, this section will have to change accordingly.
    """
    # check for empty scene configuration and set it up properly if it is empty.
    if len(scene_configuration) == 0:
        scene_configuration = build_empty_scene()

    instrument = InstrumentFactory(config=instrument_configuration, webapp=webapp)
    warnings.update(instrument.warnings)

    strategy = StrategyFactory(instrument, config=strategy_configuration, webapp=webapp)

    return calc_config, instrument, strategy, scene_configuration, background, background_level, warnings

def make_observation(calc_config, instrument, strategy, scene, background, background_level, webapp=False):
    # set up the observation and then do S/N calculation...
    obs = observation.Observation(
        scene=scene,
        instrument=instrument,
        strategy=strategy,
        background=background,
        background_level=background_level,
        webapp=webapp
    )

    # seed the random number generator
    seed = obs.get_random_seed()
    np.random.seed(seed=seed)

    # Sometimes there is more than one exposure involved so implement lists for signal and noise
    my_detector_signal_list = []
    my_detector_noise_list = []
    my_detector_saturation_list = []

    if hasattr(strategy, 'dithers'):
        dither_list = strategy.dithers
    else:
        dither_list = [{'x': 0.0, 'y': 0.0}]

    for dither in dither_list:
        o = copy.deepcopy(obs)
        o.scene.offset(dither)
        # Calculate the signal rate in the detector plane. If they're configured, need to loop through
        # configured orders to include all dispersed signal.
        if instrument.projection_type in ['slitless', 'multiorder']:
            if 'orders' in instrument.disperser_config[instrument.instrument['disperser']]:
                if instrument.instrument['filter'] in instrument.disperser_config[instrument.instrument['disperser']]['orders']:
                    orders = instrument.disperser_config[instrument.instrument['disperser']]['orders'][instrument.instrument['filter']]
                else:
                    orders = instrument.disperser_config[instrument.instrument['disperser']]['orders']
            else:
                orders = None
            if 'background_orders' in instrument.disperser_config[instrument.instrument['disperser']]:
                if instrument.instrument['filter'] in instrument.disperser_config[instrument.instrument['disperser']]['background_orders']:
                    background_orders = instrument.disperser_config[instrument.instrument['disperser']]['background_orders'][instrument.instrument['filter']]
                else:
                    background_orders = instrument.disperser_config[instrument.instrument['disperser']]['background_orders']
                background_orders = instrument.disperser_config[instrument.instrument['disperser']][
                    'background_orders']
            else:
                background_orders = []
        else:
            orders = None

        if orders is not None:
            order_signals = []
            # First, take care of the orders for which a signal is desired.
            for order in orders:
                order_signals.append(DetectorSignal(o, calc_config=calc_config, webapp=webapp, order=order))
            # Then loop through a list of orders for which we only want to add the background calculation
            if background_orders is not None:
                o_empty = copy.deepcopy(obs)
                o_empty.scene = Scene(input=build_empty_scene(), webapp=webapp)
                for order in background_orders:
                    order_signals.append(
                        DetectorSignal(o_empty, calc_config=calc_config, webapp=webapp, order=order,
                                        empty_scene=True))

            my_detector_signal = CombinedSignal(order_signals)
        else:
            my_detector_signal = DetectorSignal(o, calc_config=calc_config, webapp=webapp, order=None)

        my_detector_noise = DetectorNoise(my_detector_signal, o)

        # Every dither has a saturation map
        my_detector_saturation = my_detector_signal.get_saturation_mask()

        my_detector_signal_list.append(my_detector_signal)
        my_detector_noise_list.append(my_detector_noise)
        my_detector_saturation_list.append(my_detector_saturation)

    return my_detector_signal_list, my_detector_noise_list, my_detector_saturation_list, obs

def configure_coronagraphy(input, strategy, scene_configuration, warnings, webapp=False):
    """
    Configure the engine with the various pieces needed for a coronagraphy contrast calculation.
    """

    strategy.do_contrast = False

    # Check for user-specified dithers (the contrast calculation will add an additional 2 fictional dithers).
    if not hasattr(strategy, 'dithers') or len(strategy.dithers) != 1:
        message = "Contrast calculations currently require a single dither " \
                  "to be passed in the strategy, {} was passed".format(strategy.dithers)
        raise EngineInputError(value=message)

    psf_subtraction_configuration = input['strategy']['psf_subtraction_source']
    psf_subtraction_xy = strategy.psf_subtraction_xy
    pointing_error = strategy.pointing_error
    psf_subtraction_configuration['position']['x_offset'] = strategy.psf_subtraction_xy[0]
    psf_subtraction_configuration['position']['y_offset'] = strategy.psf_subtraction_xy[1]

    # add the psf_reference to the scene
    scene_configuration.append(psf_subtraction_configuration)
    scene = Scene(input=scene_configuration, webapp=webapp)
    if hasattr(strategy, "scene_rotation"):
        scene.sources = strategy.rotate(scene.sources)
    warnings.update(scene.warnings)
    # add an attribute specifying the PSF subtraction source
    scene.psf_subtraction_idx = len(scene_configuration) - 1

    # Add the appropriate dithers for the PSF reference star and the unocculted dither
    psf_subtraction_dither = {
        'x': -strategy.psf_subtraction_xy[0] - strategy.pointing_error[0],
        'y': -strategy.psf_subtraction_xy[1] - strategy.pointing_error[1]
    }
    unocculted_dither = {'x': strategy.unocculted_xy[0], 'y': strategy.unocculted_xy[1]}

    strategy.dithers.append(psf_subtraction_dither)
    strategy.dithers.append(unocculted_dither)
    strategy.on_target = [True, False, False]

    return scene, strategy

def calculate_sn(input, webapp=False):
    """
    This is a function to do the 'forward' exposure time calculation where, given a dict
    in engine API input format, we calculate the resulting Signal/Noise and return a Report
    on the results.

    Parameters
    ----------
    input: dict
        Engine API format dictionary containing the information required to perform the calculation.
    webapp: bool
        Keyword activating strict checking of API items

    Returns
    -------
    report.Report instance
    """
    calc_config, instrument, strategy, scene_configuration, background, background_level, warnings = setup(input, webapp=webapp)

    # strategies can have different figures of interest that need to be calculated.
    # in most cases, S/N is what is desired. However, for coronagraphy the figure of interest
    # is sometimes the contrast that can be achieved. In this case, strategy.calc_type will be
    # set to 'contrast' and we need to run calculate_contrast.
    contrast = False
    if hasattr(strategy, "calc_type"):
        if strategy.calc_type == "contrast":
            contrast = True
        else:
            msg = "Unsupported calculation type: %s" % strategy.calc_type
            raise EngineInputError(value=msg)

    if contrast:
        scene, strategy = configure_coronagraphy(input, strategy, scene_configuration, warnings, webapp=webapp)
    else:
        scene = Scene(input=scene_configuration, webapp=webapp)
        warnings.update(scene.warnings)

    my_detector_signal_list, my_detector_noise_list, my_detector_saturation_list, obs = make_observation(calc_config, instrument, 
                                                                       strategy, scene, background, background_level, webapp=webapp)

    # Use the strategy to get the extracted signal/noise products
    extracted_sn_list = strategy.extract(my_detector_signal_list, my_detector_noise_list)
    warnings.update(extracted_sn_list[0]['warnings'])


    # selection logic: If calculation_config specifies it, do that. If not, use the predefined data value
    scatter = instrument.the_detector.scatter
    if calc_config.noise['scatter'] is not None:
        scatter = calc_config.noise['scatter']
    if scatter:
        extracted_sn_list = instrument.apply_scattering(extracted_sn_list)

    if contrast:
        extracted_sn_list, warnings = calculate_contrast(my_detector_signal_list, my_detector_noise_list, extracted_sn_list, strategy, warnings, webapp=webapp)

    # #### END calculation #### #
    r = ReportFactory(input, my_detector_signal_list, my_detector_noise_list, my_detector_saturation_list,
                        extracted_sn_list, warnings)

    return r


def calculate_contrast(my_detector_signal_list, my_detector_noise_list, extracted_sn_list, strategy, warnings, webapp=False):
    """
    This is a function to calculate the resulting contrast from the coronagraphic imaging calculation. 

    While this method is meant for coronagraphic modes, it will work also for regular imaging modes.

    Parameters
    ----------
    input: dict
        Engine API format dictionary containing the information required to perform the calculation.
    psf_library : psf_library.PSFLibrary instance
        Library of PSF files (e.g. produced by webbpsf) to be used in the calculation

    Returns
    -------
    report.Report instance
    """

    # We need a regular S/N of the target source
    extracted_sn_list = strategy.extract(my_detector_signal_list, my_detector_noise_list)
    warnings.update(extracted_sn_list[0]['warnings'])

    # Use the strategy to get the extracted contrast products
    grid = my_detector_signal_list[0].grid

    aperture = strategy.aperture_size
    annulus = strategy.sky_annulus

    # Create a list of contrast separations for which to calculate the contrast
    bounds = grid.bounds()
    ncontrast = strategy.ncontrast
    contrasts = np.zeros(ncontrast)
    contrast_separations = np.linspace(0 + aperture, bounds['xmax'] - annulus[1], ncontrast)
    contrast_azimuth = np.radians(strategy.contrast_azimuth)
    contrast_xys = [(separation * np.sin(contrast_azimuth),
                     separation * np.cos(contrast_azimuth)) for separation in contrast_separations]

    # Calculate contrast at each separation
    for i, contrast_xy in enumerate(contrast_xys):
        strategy.target_xy = contrast_xy
        extracted_list = strategy.extract(my_detector_signal_list, my_detector_noise_list)
        contrasts[i] = extracted_list[0]['extracted_noise']

    # What is the flux of the unocculted star.
    # We set the do_contrast attribute to True so the unocculted dither will be used.
    strategy.do_contrast = True
    strategy.on_target = [False, False, True]
    strategy.target_xy = strategy.unocculted_xy
    extract_unocculted_list = strategy.extract(my_detector_signal_list, my_detector_noise_list)
    extracted_sn_list.append(extract_unocculted_list[0])

    # Contrast is relative to the unocculted on-axis star.
    contrasts /= extract_unocculted_list[0]['extracted_flux']
    contrast_curve = [contrast_separations, contrasts]

    # #### END calculation #### #

    # Add the contrast curve and link relevant saturation maps to the extracted_sn dict for passing
    for x in range(len(extracted_sn_list)):
        extracted_sn_list[x]['contrast_curve'] = contrast_curve

    return extracted_sn_list, warnings

def calculate_exposure_time(input, webapp=False, **kwargs):
    """
    This is a function to do the 'reverse' exposure time calculation where given a desired
    signal-to-noise ratio we calculate the optimal exposureSpecification and return a Report
    on the results.

    Parameters
    ----------
    input: dict
        Dictionary containing the information required to perform the calculation.
    webapp: bool
        Keyword activating strict checking of API items

    Returns
    -------
    report.Report instance
    """
    calc_config, instrument, strategy, scene_configuration, background, background_level, warnings = setup(input, webapp=webapp)
    # do a forward calculation with time=1
    # this sets the detector in instrument AND returns a new input dict
    instrument.the_detector.set_time(1.0)

    # strategies can have different figures of interest that need to be calculated.
    # in most cases, S/N is what is desired. However, for coronagraphy the figure of interest
    # is sometimes the contrast that can be achieved. In this case, strategy.calc_type will be
    # set to 'contrast' and we need to run calculate_contrast.
    contrast = False
    if hasattr(strategy, "calc_type"):
        if strategy.calc_type == "contrast":
            contrast = True
        else:
            msg = "Unsupported calculation type: %s" % strategy.calc_type
            raise EngineInputError(value=msg)

    if contrast:
        scene, strategy = configure_coronagraphy(input, strategy, scene_configuration, warnings, webapp=webapp)
    else:
        scene = Scene(input=scene_configuration, webapp=webapp)
        warnings.update(scene.warnings)

    my_detector_signal_list, my_detector_noise_list, my_detector_saturation_list, obs = make_observation(calc_config, instrument, 
                                                                       strategy, scene, background, background_level, webapp=webapp)

    # Use the strategy to get the extracted signal/noise products
    extracted_sn_list = strategy.extract(my_detector_signal_list, my_detector_noise_list)
    warnings.update(extracted_sn_list[0]['warnings'])

    new_time = my_detector_signal_list[0].the_detector.scale_exposure_time(extracted_sn_list[0], input)
    obs.instrument.the_detector.set_time(new_time)

    for signal in range(len(my_detector_signal_list)):
        # reset the time(s)
        my_detector_signal_list[signal].the_detector.set_time(new_time)
        my_detector_signal_list[signal].instrument.the_detector.set_time(new_time)
        
        # retrigger calculations using the new time
        my_detector_signal_list[signal].trigger_saturation_calcs()
        my_detector_noise_list[signal] = DetectorNoise(my_detector_signal_list[signal], obs)
        my_detector_saturation_list[signal] = my_detector_signal_list[signal].get_saturation_mask()


    # Use the strategy to get the extracted signal/noise products
    extracted_sn_list = strategy.extract(my_detector_signal_list, my_detector_noise_list)

    # selection logic: If calculation_config specifies it, do that. If not, use the predefined data value
    scatter = instrument.the_detector.scatter
    if calc_config.noise['scatter'] is not None:
        scatter = calc_config.noise['scatter']
    if scatter:
        extracted_sn_list = instrument.apply_scattering(extracted_sn_list)


    if contrast:
        extracted_sn_list, warnings = calculate_contrast(my_detector_signal_list, my_detector_noise_list, extracted_sn_list, strategy, warnings, webapp=webapp)

    # #### END calculation #### #
    r = ReportFactory(input, my_detector_signal_list, my_detector_noise_list, my_detector_saturation_list,
                        extracted_sn_list, warnings)

    return r
