"""This module defines functions for converting BIT queries into pyetc
engine inputs.
"""
import os, pprint, re, sys
import six

from   pandeia.engine.helpers.bit.pyetc_util import log, read_dict
import pandeia.engine.helpers.bit.pyetc_form_defaults as pyetc_defaults

# ========================================================================

def query_to_web_inputs(query_str):
    """Translate a single query string into a dictionary of partial
    JETC web inputs and a list of requested JETC result variables.

    query_str --> parital_web_inputs, requested_jetc_results
    """
    query_dict = query_to_dict(query_str)
    bit_inputs, requests = partition_query_dict(query_dict)
    web_inputs = expand_bit_inputs(bit_inputs)
    return web_inputs, requests

# ========================================================================

def query_to_dict(query):
    """Converts a web query GET parameter string into a Python dict mapping
    strings to strings.

    >>> query_to_dict("a=1&b=2&c=3")
    {'a': '1', 'c': '3', 'b': '2'}
    """
    keyvals = query.strip().split("&")
    kvpairs = []
    for assignment in keyvals:
        kvt = tuple([x.strip() for x in assignment.split("=")])
        if len(kvt) == 1:
            kvt = kvt + ("",)
        kvpairs.append(tuple(kvt))
    return dict(sorted(kvpairs))

def partition_query_dict(query_dict):
    """query_dict -->  (bit_inputs,  requests)"""
    bit_inputs = dict(query_dict)
    requests = []
    request_count = -1
    bit_inputs_keys = sorted(bit_inputs.keys())
    for key in bit_inputs_keys:
        if key.startswith("request"):
            if key == "requestCount":
                request_count = int(bit_inputs.pop(key))
            else:
                requests.append(bit_inputs.pop(key))
    assert len(requests) == request_count, "wrong number of requests"
    return bit_inputs, requests

# ========================================================================

def expand_bit_inputs(bit_inputs, config=None):
    """Expand `bit_inputs` into pyetc web inputs using `config` spec.   If
    `config` is None,  get the spec based on the "config" query variable.

    Return pyetc web inputs.
    """
    result = dict()
    if config is None:
        config = get_bit_config(bit_inputs["config"])
    for key, val in config.items():
        if isinstance(val, (six.string_types, six.integer_types, float, bool)):
            result[key] = str(val)
        elif isinstance(val, dict) and key in bit_inputs:
            case = bit_inputs[key]
            for cases, what_to_do in val.items():
                if case in cases:
                    result.update(expand_bit_inputs(bit_inputs, what_to_do))
                    break
            else:
                raise ValueError("No handler found for " + key + " = " +
                                 repr(bit_inputs[key]))
        else:
            raise RuntimeError("Specification error at " + repr((key,val)))
    expanded = subst_vars(bit_inputs, result)
    return expanded

HERE = os.path.dirname(__file__) or "./"

def load_bit_config(config):
    """Load a single BIT configuration file combined with global settings.
    """
    try:
        path = HERE + "/config/global.dict"
        conf = read_dict(path)
        path = HERE + "/config/" + config.replace("/","_").lower() + ".dict"
        conf.update(read_dict(path))
    except IOError:
        raise ValueError("Unknown query config = " + repr(config))
    return conf

CONFIG_CACHE = {}
def get_bit_config(config):
    """Get a query expansion spec for a particular `config` setting,  looking
    first in the CONFIG_CACHE and loading the `config` only if not found.

    >>> _ = get_bit_config("acs/wfc")
    >>> _ = get_bit_config("acs/sbc")
    >>> _ = get_bit_config("cos/nuv")
    >>> _ = get_bit_config("cos/fuv")
    >>> _ = get_bit_config("wfc3/ir")
    >>> _ = get_bit_config("wfc3/uvis")
    >>> _ = get_bit_config("stis/ccd")
    >>> _ = get_bit_config("stis/nuv-mama")
    >>> _ = get_bit_config("stis/fuv-mama")
    """
    if config not in CONFIG_CACHE:
        CONFIG_CACHE[config] = load_bit_config(config)
    return CONFIG_CACHE[config]

LC_PERC_VAR = re.compile("%%([A-Za-z0-9_]+)%%")
PERC_VAR = re.compile("%([A-Za-z0-9_]+)%")

def subst_var(bit_inputs, string):
    """Replace %var% in `string` with bit_inputs[var].
    Replace %%var%% in `string` with bit_inputs[var].lower().
    Otherwise return `string`.
    """
    match = LC_PERC_VAR.search(string)
    if match:
        return LC_PERC_VAR.sub(bit_inputs[match.group(1)].lower(),
                           string)
    match = PERC_VAR.search(string)
    if match:
        return PERC_VAR.sub(bit_inputs[match.group(1)],
                            string)
    return string

def subst_vars(bit_inputs, config):
    """Substitute values from `bit_inputs` into all strings in BIT `config`,
    whether they be keys, values, simple strings, elements of lists, etc.
    """
    if isinstance(config, six.string_types):
        result = subst_var(bit_inputs, config)
    elif isinstance(config, (list, tuple)):
        result = []
        for elem in config:
            result.append(subst_vars(bit_inputs, elem))
        result = tuple(result)
    elif isinstance(config, dict):
        result = {}
        for key, value in config.items():
            result[subst_vars(bit_inputs, key)] = subst_vars(bit_inputs, value)
    elif isinstance(config, (six.integer_types, float, bool)):
        result = str(config)
    else:
        raise TypeError("Don't know how to handle " + repr(config))
    return result

# ========================================================================

def convert_web_inputs(web_inputs):
    """Converts partial pyetc web inputs into fully specified pyetc
    web inputs.
    """
    instrument, sciencemode = get_migration_instr_sciencemode(web_inputs)

    web_inputs["instrument"] = web_inputs["instrument"].lower()
    web_inputs["science_mode"] = web_inputs["science_mode"].lower()

    # Get defaults for this instrument and sciencemode from the scraped
    # values kept by the spider tools.
    converted = pyetc_defaults.get_defaults(instrument, sciencemode)
    converted.update(web_inputs)

    # if instrument == "cos":
    #     hack_cos_inputs(converted)

    return converted

def hack_cos_inputs(jraw_inputs):
    """Do in-place hacks to fix wavelength issues for COS centralWavelength
    and obswave.
    """
    # This j2p hack may belong in the migration tools which are currently
    # supplying the obsolete JETC default of 1230.
    g140l = jraw_inputs.get("G140L_CentralWavelength", None)
    if g140l == "1230":
        jraw_inputs["G140L_CentralWavelength"] = "1280"

    # Ordinarily,  obswave is assigned the same value as centralWavelength from the query.
    # Sometimes it doesn't work so we hack just the required points until the full set of
    # ranges is implemented somewhere in Python.
    # RangeError("Wavelength 1105.0 not found in available cos segment ranges {(1120.0, 2250.0): 'a'}",)
    # RangeError("Wavelength 1280.0 not found in available cos segment ranges {(1281.0, 2385.0): 'a', (700.0, 1195.5): 'b'}",)
    # RangeError("Wavelength 1300.0 not found in available cos segment ranges {(1146.0, 1287.0): 'b', (1301.0, 1442.0): 'a'}",)
    # RangeError("Wavelength 1600.0 not found in available cos segment ranges {(1410.0, 1581.0): 'b', (1601.0, 1772.0): 'a'}",)
    hacked_obswaves = {
        "1105" : "1120",
        "1280" : "1281",
        "1300" : "1301",
        "1600" : "1601",
    }
    if "obswave" in jraw_inputs:
        obswave = jraw_inputs.get("obswave")
        jraw_inputs["obswave"] = hacked_obswaves.get(obswave, obswave)

def get_migration_instr_sciencemode(web_inputs):
    """Determine the instrument and science_mode with the spelling used by
    the spider tools.
    """
    instrument = web_inputs["instrument"].lower()
    sciencemode = {
        "Imaging" : "imag",
        "Spectroscopic" : "spec",
        "Rampfilter" : "ramp"
    }.get(web_inputs["science_mode"])
    return instrument, sciencemode

# ========================================================================

# Function to help analyze queries to determine parameter space.

def get_query_values(filename):
    """Utility to return the possible values for each variable defined in
    some query in `filename`.   For development analysis.

    queries(filename) -->  { var : [values...], ... }
    """
    values = {}
    for query in queries(filename):
        query_dict = query_to_dict(query)
        for key, value in query_dict.items():
            if key not in values:
                values[key] = [value]
            elif value not in values[key]:
                values[key].append(value)
    return values

# ========================================================================

# The main translation function called by BIT computations.

def process_query(query_str, info=True):
    """Translate and dump a single query.

    Returns {web_inputs}, [requested_outputs]
    """
    if info:
        log.info("BIT Query:", repr(query_str))
    web_inputs, requests = query_to_web_inputs(query_str)
    complete_web_inputs = convert_web_inputs(web_inputs)
    if log.VERBOSE_FLAG:
        log.verbose("\nComplete pyetc web inputs:\n",
                    pprint.pformat(complete_web_inputs))
    return complete_web_inputs, requests


# ========================================================================

# A regression test function to determine that engine inputs haven't changed.

def process_file(fname):
    """Translate an entire file containing one query per line.
    Outputs to a like-named .qweb file:  { query : { web_inputs } }
    """
    ofname = os.path.splitext(fname)[0]+".qweb"
    ofile = open(ofname, "w+")
    ofile.write("{\n")
    for query_str in queries(fname):
        try:
            web_inputs, requests = process_query(query_str, info=True)
        except Exception:
            log.error("Exception in " + repr(query_str))
        else:
            web_inputs["bit_requests"] = requests
            ofile.write(repr(query_str) + " : " + repr(web_inputs) + ",\n\n")
    ofile.write("}\n")

def queries(filename):
    """Generator returning the query strings defined in `filename`,
    removing comments.
    """
    for line in open(filename):
        if not line.strip().startswith("//"):
            yield line.strip()
