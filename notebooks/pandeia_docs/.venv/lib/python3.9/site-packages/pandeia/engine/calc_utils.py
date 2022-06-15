# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, print_function, absolute_import

import os
import pkg_resources
import numpy as np

from . import config
from .io_utils import read_json
from pandeia.engine.source import Source, Line
from pandeia.engine.custom_exceptions import EngineInputError
from pandeia.engine.normalization import NormalizationFactory
from pandeia.engine.sed import SEDFactory

default_refdata_directory = config.default_refdata_directory

# Python 2 triggers IOError instead of FileNotFoundError.
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = OSError


def strip_meta(d):
    """
    recursively strip any 'meta' tags from a dict, d. this function modifies the dict in place.
    """
    if isinstance(d, dict):
        if 'meta' in d:
            d.pop('meta')
        for k in d:
            strip_meta(d[k])


def get_telescope_config(telescope):
    """
    Get the telescope configuration to tell us instruments and modes.
    Must have the pandeia_refdata env variable set and pointing to the reference data

    Parameters
    ----------
    telescope: string
        Telescope identifier string. Currently supports jwst, roman, and hst.

    Returns
    -------
    config: dict
        Telescope configuration data
    """
    config_file = os.path.join(default_refdata_directory, telescope, 'telescope', 'config.json')
    config = read_json(config_file, raise_except=True)
    strip_meta(config)
    return config


def get_instrument_config(telescope, instrument):
    """
    Get the instrument configuration.
    Must have the pandeia_refdata env variable set and pointing to the reference data

    Parameters
    ----------
    telescope: string
        Telescope identifier string. Currently supports jwst, roman, and hst.
    instrument: string
        Instrument identifier string.

    Returns
    -------
    config: dict
        Instrument configuration data
    """
    config_file = os.path.join(default_refdata_directory, telescope, instrument, 'config.json')
    config = read_json(config_file, raise_except=True)
    strip_meta(config)
    return config


def build_default_source(geometry="point"):
    """
    Build and return a default Source in engine API dict format

    Returns
    -------
    src: dict
        Source configuration data in engine API format
    """
    src = Source(shape={"geometry": geometry}).as_dict()
    src['spectrum']['normalization'] = NormalizationFactory(config=src['spectrum']['normalization']).as_dict()
    src['spectrum']['sed'] = SEDFactory(config=src['spectrum']['sed']).as_dict()
    return src


def build_default_line():
    """
    Build and return a default Line in engine API dict format

    Returns
    -------
    line: dict
        Line configuration data in engine API format
    """
    line = Line().as_dict()
    return line


def build_default_scene():
    """
    Build and return a default scene which consists of a single default Source

    Returns
    -------
    scene: list
        Single element list containing a source as built by build_default_source()
    """
    scene = [build_default_source()]
    return scene


def build_empty_scene():
    """
    Build and return an empty scene. Because of the way ConvolvedSceneCube's and AstroSpectrum's are created,
    there must be some concept of a spectrum and thus a source contained within a scene. The (admittedly hacky)
    workaround is to build a default scene with a flat-spectrum point source and set its flux to 0. This will provide
    the framework necessary to build things from this scene using the existing code, but not add any actual flux. The
    primary motivation for this is the BNC which needs to build a ConvolvedSceneCube that only contains background signal.

    Returns
    -------
    scene: list
        Single element list containing a zero-flux source
    """
    scene = build_default_scene()
    scene[0]['spectrum']['normalization']['norm_flux'] = 0.0
    return scene


def build_default_calc(telescope, instrument, mode, **kwargs):
    """
    Build a default calculation given an telescope, instrument, and mode

    Parameters
    ----------
    telescope: string
        Telescope identifier string. Currently supports jwst, roman, and hst.
    instrument: string
        Instrument identifier string.
    mode: string
        Mode identifier string.

    Returns
    -------
    calc: dict
        Engine output API compliant dict defining a calculation
    """
    calc = dict()
    ex_args = dict(**kwargs)

    ins_config = get_instrument_config(telescope=telescope, instrument=instrument)
    if mode is None or mode not in ins_config['modes']:
        mode = ins_config['defaults']['mode']

    # the instrument configuration now contains instrument and strategy values
    calc['configuration'] = ins_config['defaults'][mode]['instrument_args']

    calc['scene'] = build_default_scene()

    bg_defaults_file = "defaults/background_defaults.json"
    bg_defaults = read_json(pkg_resources.resource_filename("pandeia.engine", bg_defaults_file))

    calc['background'] = bg_defaults[telescope]['bg_location']
    calc['background_level'] = bg_defaults[telescope]['bg_level']

    calc_defaults_file = "defaults/calculationconfig_defaults.json"
    calc['calculation'] = read_json(pkg_resources.resource_filename("pandeia.engine", calc_defaults_file))

    st_defaults = ins_config['strategy_config'][mode]
    calc['strategy'] = dict()

    if 'method' in ex_args:
        calc['strategy']['method'] = ex_args['method']
    else:
        calc['strategy']['method'] = st_defaults['method']

    # Only use the constructed strategy with the default if this IS the default strategy.
    if ins_config['defaults'][mode]['strategy_args']['method'] == calc['strategy']['method']:
        calc['strategy'] = ins_config['defaults'][mode]['strategy_args']
    else:
        # otherwise, construct the entire strategy
        method_defaults_file = "%s.json" % calc['strategy']['method']

        method_defaults = read_json(os.path.join(default_refdata_directory, "strategy", method_defaults_file))
        for k in method_defaults:
            if 'permitted' not in k:
                calc['strategy'][k] = method_defaults[k]

        # sometimes aperture_key is a dict because some strategies will use it, but others won't.
        # if it's not used, it will be set to None
        if 'aperture_key' in st_defaults:
            if isinstance(st_defaults['aperture_key'], dict):
                ap_key = st_defaults['aperture_key'][calc['strategy']['method']]
            else:
                ap_key = st_defaults['aperture_key']
            if ap_key is not None:
                key = calc['configuration']['instrument'][ap_key]
                calc['strategy']['aperture_size'] = st_defaults['aperture_sizes'][key]
        else:
            ap_key = None
            if 'aperture_size' in st_defaults:
                calc['strategy']['aperture_size'] = st_defaults['aperture_size']

        # in most cases sky_key will be the same as ap_key.  however, in the case of IFUNodApPhot, aperture
        # is used but NOT sky_annulus.  this is here to deal with that scenario, but may come up in future
        # scenarios as well.  if sky_annulus is not to be used, then sky_key will be set to None.

        #  NOTE this should get refactored to be more generalized.  use of sky_key for coronagraphy
        #  parameters is hackish...
        if "sky_key" in st_defaults:
            sky_key = st_defaults['sky_key'][calc['strategy']['method']]
        else:
            sky_key = ap_key

        # make sure sky_key is set and sky_annulus_sizes are actually defined
        params = {
            "sky_annulus_sizes": "sky_annulus",
            "unocculted_xys": "unocculted_xy",
            "contrast_azimuths": "contrast_azimuth",
            "contrast_separations": "contrast_separation"
        }
        for p in params:
            if sky_key is not None and p in st_defaults:
                key = calc['configuration']['instrument'][sky_key]
                calc['strategy'][params[p]] = st_defaults[p][key]

    return calc


def calcs_to_wb(calc_list, name="Generated workbook", note="", desc="", wb_id=1, proposal_id="", proposal_state=""):
    """
    Take a list of calculations and make a workbook out of them.  This routine doesn't do any redundancy checking
    between scenes and sources so each calculation has its own scene with its own distinct sources.

    Parameters
    ----------
    calc_list: list of dicts (engine input API)
        List of dicts in engine input API compatible format

    Returns
    -------
    wb: dict
        Workbook data in pandeia workbook format
    """
    wb = {}
    scenes = {}
    sources = {}
    source_i = 1
    calcs = {}
    calc_i = 1
    for c in calc_list:
        # set up the scene for this calculation
        scenes[str(calc_i)] = {}
        scenes[str(calc_i)]["deleted"] = 0
        scenes[str(calc_i)]["characteristics"] = "{}"
        scenes[str(calc_i)]["name"] = "Scene #%d" % calc_i
        scenes[str(calc_i)]["desc"] = "Scene for calculation #%d" % calc_i
        scenes[str(calc_i)]["scene_elements"] = {}
        for src in c['scene']:
            # go through the scene and add the sources there to the workbook's list of sources
            scenes[str(calc_i)]["scene_elements"][str(source_i)] = src['position']
            scenes[str(calc_i)]["scene_elements"][str(source_i)]['name'] = "%d Calculation %s" % (source_i, calc_i)
            scenes[str(calc_i)]["scene_elements"][str(source_i)]['deleted'] = 0
            scenes[str(calc_i)]["scene_elements"][str(source_i)]['source_id'] = source_i
            scenes[str(calc_i)]["scene_elements"][str(source_i)]['source_arg_blob'] = []
            sources[str(source_i)] = {}
            sources[str(source_i)]['deleted'] = 0
            sources[str(source_i)]['name'] = "%f %s %s source" % (
                src['spectrum']['normalization']['norm_flux'],
                src['spectrum']['normalization']['norm_fluxunit'],
                src['shape']['geometry']
            )
            sources[str(source_i)]['characteristics'] = {}
            sources[str(source_i)]['characteristics']['shape'] = src['shape']
            sources[str(source_i)]['characteristics']['spectrum'] = src['spectrum']
            source_i += 1
        # set up the calculations
        calcs[str(calc_i)] = {}
        calcs[str(calc_i)]["strategy_args"] = c['strategy']
        calcs[str(calc_i)]["name"] = "Calculation #%d" % calc_i
        calcs[str(calc_i)]["deleted"] = 0
        calcs[str(calc_i)]["apt"] = 0
        calcs[str(calc_i)]["strategy"] = c['strategy']['method']
        calcs[str(calc_i)]["client_data"] = {
            "cause": "generated_wb"
        }
        calcs[str(calc_i)]["background"] = {"bg_type": c['background'], "ra": 0.0, "dec": 0.0}
        calcs[str(calc_i)]["scene_id"] = calc_i
        calcs[str(calc_i)]["camera_config"] = c['configuration']
        calcs[str(calc_i)]["instrument"] = calcs[str(calc_i)]["camera_config"]["instrument"].pop("instrument")
        calcs[str(calc_i)]["insmode"] = calcs[str(calc_i)]["camera_config"]["instrument"].pop("mode")
        calc_i += 1
    # now populate the workbook
    wb["test_mode"] = "1"
    wb["created"] = 0
    wb["deleted"] = 0
    wb["id"] = wb_id
    wb["name"] = name
    wb["desc"] = desc
    wb["note"] = note
    wb["proposal_id"] = proposal_id
    wb["proposal_state"] = proposal_state
    wb["calculations"] = calcs
    wb["scenes"] = scenes
    wb["sources"] = sources
    return wb
