""" SEE NOTES IN pyetc repo  """

import pprint as pp
from xml.sax import saxutils

from pandeia.engine.helpers.bit.pyetc_util import log
from pandeia.engine.helpers.bit import input_conversion
from pandeia.engine.helpers.bit.perform_calculation import process_form_input

# ========================================================================

OUTPUT_MAP = {
    "BrightestPixelRate" : "brightest_pixel_rate",
    "ImageRate" : "detector_total_rate",
    "BackgroundRate" : "extracted_background_rate",

     'SegmentedCountRate.A': 'segment_a_total_rate',
     'SegmentedCountRate.B': 'segment_b_total_rate',
     'SegmentedCountRate.C': 'stripe_c_total_rate',

     # Forecasted mappings from Meg based on unimplemented pyetc functionality
     'SegmentedBrightPixel.A.value' : \
        lambda v, r: r["segmented_brightest_pixel_rate"]["a"][0],
     'SegmentedBrightPixel.B.value' : \
        lambda v, r: r["segmented_brightest_pixel_rate"]["b"][0],
}

def _has_segment_b(ei):
    """ Returns True if this is a COS calculation for which
    we would expect results for Segment B (not all do).
    Checks over the engine inputs (`ei`) dict passed in.
    """
    # sanity checks
    if ei['instrument'] != 'cos': raise ValueError('_has_segment_b called for non-COS calculation')
    for item in ('detector', 'grating', 'central_wavelength'):
        if item not in ei:
            raise ValueError(item+' not in engine inputs')

    segaltnames_key = (ei['grating'], ei['central_wavelength'])
    # check the COS instrument data config dict.  if segaltnames_key is in
    # it, then that will tell us if we should expect segment B.  if not, assume
    # it has a segment B
    detector = ei['detector']
    import pyetc.instruments.hst.cos as pihc # this is cached after 1st call
    cos_det_cfg = pihc.data.dconfig[detector]
    if segaltnames_key in cos_det_cfg['segaltnames']:
        segnames = cos_det_cfg['segaltnames'][segaltnames_key]
        return 'b' in segnames or 'B' in segnames
    else:
        return True

def get_and_map_requested_results(query, full_results, requests):
    """Map the variable names in `requests` from JETC nomenclature
    onto pyetc nomenclature and return them.  Return unavailable
    request variables as "undefined". The `query` var is passed in
    only to aid in debugging.

    dict(full_pyetc_results), list(requests) -->  dict(requested_jetc_results)
    """
    partial_results = {}
    for var in requests:
        mapped_var = None
        try:
            mapped_var = OUTPUT_MAP[var]
            # obsmode = full_results["engine_inputs"]["obsmode"]
            # grating = full_results["engine_inputs"].get("grating","none")
            if isinstance(mapped_var, str):
                partial_results[var] = full_results[mapped_var]
            else:
                partial_results[var] = mapped_var(var, full_results)
        except KeyError:
            partial_results[var] = "undefined"
            if mapped_var == None:
                log.error("Missing requested var from OUTPUT_MAP:", var)
            else:
                # The mapped_var is missing from full_results, so error, but if
                # this is segment B for COS, first see if we should even have it
                if full_results['engine_inputs']['instrument'] == 'cos' and \
                   var.startswith('Segmented') and '.B' in var and \
                   not _has_segment_b(full_results['engine_inputs']):
                    pass # this is OK - there is no segment B data
                else:
                    log.error("Missing requested mapped var from full results: "+str(mapped_var))
                    log.error("... query was:\n"+query+"\n... full_results is:\n"+str(full_results)+"\n")
    return partial_results

# ========================================================================

def get_partial_results(query, requests, full_results, message):
    """Return an XML <Response> associated with the result of a query.

    Pull out JETC named `requests` from pyetc:JETC mapped `full_results`.

    Include `query` in as the <query> element.

    Set <validRequest> to True IFF message is None.  Otherwise set
    <validRequest> to False and add <message>.

    Return a string.
    """
    partial_results = {}
    if message is not None:
        partial_results["message"] = message
        partial_results["validRequest"] = "false"
    else:
        requested = get_and_map_requested_results(query, full_results, requests)
        partial_results.update(requested)
        partial_results["validRequest"] = "true"
    partial_results["query"] = query
    xml = dict_to_xml("Response", partial_results)
    return str(xml) + "\n"

def dict_to_xml(name, dict_):
    """Serialize `dict_` as an XML element named `name`.

    >>> dict_to_xml("top_level", {
    ...        "this" : {
    ...            "nested":1,
    ...            "other":"str"
    ...        },
    ...        "that": 42
    ...    })
    '<top_level><this><other>str</other><nested>1</nested></this><that>42</that></top_level>'
    """
    return _make_element(name, _value_to_xml(dict_))

def _make_element(name, innards):
    """Wrap string `innards` in an XML begin/end tag."""
    return "<" + name + ">" + innards + "</" + name + ">"

def _value_to_xml(item):
    """Convert value `item` into XML,  recursively converting dictionaries
    into nested XML and XML-escaping the str() of non-dictionaries.
    """
    if isinstance(item, dict):
        xml = []
        for key, val in list(item.items()):
            xml.append(_make_element(key, _value_to_xml(val)))
        return "".join(xml)
    else:
        return saxutils.escape(str(item))

# ========================================================================

def compute_bright_object(web_inputs):
    """Given a set of pyetc web_inputs,  return pyetc .r dict with
    extra "engine_inputs" field.
    """
    engine_inputs = web_to_engine_inputs(web_inputs)

    results = enginecompute.compute(engine_inputs.copy(),
                                    safe_exceptions=False)
    results["engine_inputs"] = engine_inputs
    return results


def web_to_engine_inputs(web_inputs):
    """Convert complete pyetc web inputs into engine inputs.
    """
    form_dict = input_conversion.convert_dictionary(web_inputs)[0]
    instr = form_dict["instrument"]
    scimode = form_dict["science_mode"]

    engine_inputs = process_form_input(instr, scimode, form_dict, "./")[0]

    if log.VERBOSE_FLAG:
        log.verbose("\npyetc engine inputs:\n", pp.pformat(engine_inputs), "\n")

    return engine_inputs
