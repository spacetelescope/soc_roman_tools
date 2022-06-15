#! /usr/bin/env python
""" This module is the main script for BIT batch mode processing. This
is an edited version of the original copied from the pyetc codebase.

The plan is it takes an input file with one query per line and produces
an output file with one result per line.

It WILL BE structured to
facilitate parallel processing but for now can only be run serially.
"""
import io, optparse, os, pprint, sys, json
import signal
from multiprocessing import Pool
from datetime import datetime
import pysynphot

from pandeia.engine.helpers.bit.pyetc_util import log
from pandeia.engine.helpers.bit import translate
from pandeia.engine.helpers.bit import compute
from pandeia.engine.helpers.peng import peng3jeng
from pandeia.engine.perform_calculation import perform_calculation

#from pyetc.etc_web.etc import versions

# ========================================================================

# process_query is the main parallelizable function.  It
# computes an ETC for one query at a time.  Since queries are
# independent, many can be run simultaneously.  Global references need
# to be handled specially for IPython's parallel processing: imports
# and cache.  Likewise exception handling is designed to trap all
# exceptions and funnel them back as either a response <message> or an
# exception string which is re-raised by the top-level loop, i.e. the
# serial controller.

RESULT_2_XML = {
    'brightest_pixel_rate': 'BrightestPixelRate',
    'detector_total_rate': 'ImageRate'
}

def make_element(name, value):
    return '<' + name + '>' + str(value) + '</' + name + '>'

def get_version_string():
    #TODO: pandeia doesn't have a __version__ attribute, so it's hardcoded here.  Probably a better way to handle this.
    versions = '<VersionInfo>\n'
    versions += '    ' + make_element('Pandeia', 'pandeia_1.7') + '\n'
    versions += '    ' + make_element('Bit', 'pandeia_1.7') + '\n'
    versions += '    ' + make_element('Pysynphot', pysynphot.__version__) + '\n'
    versions += '    ' + make_element('Cdbs', os.environ['PYSYN_CDBS']) + '\n'
    versions += '    ' + make_element('DateProcessed', datetime.now().strftime('%Y-%m-%d %H:%M:%S')) + '\n'
    versions += '</VersionInfo>'

    return versions

def create_xml_response(query, output_dict, detectors):
    # Grab the required values from the output dictionary.
    # TODO: This works for STIS, but for instruments with more than 1 detector we will also need to read 'b', 'c', etc
    bpr_name = 'brightest_pixel_rate'
    ir_name = 'detector_total_rate'
    bpr_xml_name = RESULT_2_XML[bpr_name]
    ir_xml_name = RESULT_2_XML[ir_name]

    found_detector = False
    for detector in detectors:
        if detector in output_dict['bop'][bpr_name]:
            bpr = output_dict['bop'][bpr_name][detector][0]
            ir = output_dict['bop'][ir_name][detector]
            found_detector = True

    if not found_detector:
        raise Exception("Calculation includes detector not found in detectors list.  May need to update the detectors list in bit.py")

    #bpr_wavelength = output_dict['bop'][bpr_name]['a'][1]
    #ir_wavelength = output_dict['bop'][ir_name]['a'][1]
    vr = "true"

    xml = '<Response>'
    xml += make_element(bpr_xml_name, bpr)
    xml += make_element(ir_xml_name, ir)
    xml += make_element('validRequest', vr)

    # Escape the ampersands and strip the new line character from the end of the query.
    #  The new line was causing the failure in ingestion.
    xml += make_element('query', query.replace('&', '&amp;').strip())
    xml += '</Response>'

    return xml

def get_detectors_list():
    data_dir = os.environ['pandeia_refdata']
    detectors = []

    # Open the HST telescope configuration file to read all of the instruments.
    with open(data_dir + '/hst/telescope/config.json') as tele_file:
        tele_dict = json.load(tele_file)

        # Loop through the instruments, and get all of the detectors.
        for instrument in tele_dict['instruments']:
            with open(data_dir + '/hst/%s/config.json' % instrument) as inst_file:
                inst_dict = json.load(inst_file)
                detectors.extend(inst_dict['detectors'])

    return detectors

def _process_query(query_str):
    """process_query() translates a query_str into an ETC,  runs it,
    and returns an XML response string containing requested variables.
    This is the main function designed for parallel processing,  hence
    it uses no global variables,  not even modules.
    """
    # Don't catch <control>-C in child processes
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    # First translate and run the ETC, returning exceptions as
    # <message> text.
    try:
        web_inputs, requests = translate.process_query(query_str, info=False)
        full_results = compute.compute_bright_object(web_inputs)
        message = None
    except Exception as excp:
        full_results = None
        requests = []
        message = repr(excp)

    # Next extract the requested variables from the results and format
    # XML.  Note that any exception in this phase prevents the
    # creation of the XML and <message> element used to report
    # exceptions above.  So these are actually Response-less XML-less
    # results.  Generally there should be no exceptions here since
    # even undefined request variables are already trapped and
    # handled.
    try:
        requested_xml = compute.get_partial_results(
            query_str, requests, full_results, message)
        return True, query_str, requested_xml
    except Exception as excp:   # Didn't get XML
        return False, query_str, str(excp)

def process_query(query_str):
    """Capture subprocess output to prevent scrambling during
    multiprocessing.
    """
    output = io.StringIO()
    sys.stdout, sys.stderr = output, output
    result = _process_query(query_str)
    sys.stdout, sys.stderr = sys.__stdout__, sys.__stderr__
    return result, output.getvalue().strip()

# ========================================================================

def version_string(subsystem):
    """Return a string identifying the version of `subsystem` being used."""
    if subsystem == "cdbs":
        return os.environ["PYSYN_CDBS"]
    else:
        retval = ''
        vers = versions.get_module_versions()
        if subsystem == 'bit':
            subsystem = 'pyetc'
        # get basic VCS revision
        retval = vers[subsystem]["rev"]
        # handle pyetc.* subsystems which no longer have their own version
        if subsystem.startswith('pyetc.') and retval == '':
            retval = vers['pyetc']["rev"]
        # handle non-pyetc subsystems which only have ver info in the 'str' field
        if retval == '' and vers[subsystem]['str'] not in ('', 'unknown'):
            return vers[subsystem]['str']
        # prepend version name if any
        if subsystem.startswith('pyetc') and vers[subsystem]['file'] not in ('', 'unknown'):
            vername = [pth for pth in vers[subsystem]['file'].split('/') if pth.startswith('pyetc_')][0]
            retval = vername+' '+retval
        # finally
        return retval

def get_and_log_version(subsystem):
    """Get a software version, return it, and report it in the log."""
    vers = version_string(subsystem)
    log.info("Version", subsystem, vers)
    return vers

def process_queries(inname, outname, n_compute_nodes=1):
    """process_queries() defines the top level loop over all the queries
    found in file `inname`,  writing XML results to file `outname`.
    When n_compute_nodes > 1,  process_queries is designed to parallel
    process small numbers of queries via the `mapper` function.
    """
    log.info("Processing with", n_compute_nodes, "compute nodes.")

    engine_vers = get_and_log_version("pyetc.engine")
    instruments_vers = get_and_log_version("pyetc.instruments")
    web_vers = get_and_log_version("pyetc.etc_web")
    bit_vers = get_and_log_version("bit")
    pysynphot_vers =get_and_log_version("pysynphot")
    cdbs_vers = get_and_log_version("cdbs")

    pool = Pool(n_compute_nodes)

    simple_file = open(outname,"w+")
    simple_file.write("<Responses>\n")

    simple_file.write("<VersionInfo>\n")
    simple_file.write("    <PyetcEngine>" + engine_vers + "</PyetcEngine>\n")
    simple_file.write("    <PyetcInstruments>" + instruments_vers + "</PyetcInstruments>\n")
    simple_file.write("    <PyetcWeb>" + web_vers + "</PyetcWeb>\n")
    simple_file.write("    <Bit>" + bit_vers + "</Bit>\n")
    simple_file.write("    <Pysynphot>" + pysynphot_vers + "</Pysynphot>\n")
    simple_file.write("    <Cdbs>" + cdbs_vers + "</Cdbs>\n")
    simple_file.write("    <DateProcessed>" +
        str(datetime.datetime.now()).split(".")[0] + "</DateProcessed>\n")
    simple_file.write("</VersionInfo>\n")

    # use 1 as a sentinel to tell us not to use any multiprocessing code
    parallel_process = n_compute_nodes > 1

    # chunks() returns n_compute_nodes-length lists generated by generator
    # processing by chunks helps generate incremental output.
    query_generator = translate.queries(inname)
    if not parallel_process:
        for query in query_generator:
            log.info("BIT Query:", repr(query))
            (xmlok, query, requested_xml), output = process_query(query)
            log.write(output, eol="")
            if xmlok:
                simple_file.write(requested_xml)
            if not xmlok or "<message>" in requested_xml:
                log.error("Exception:", requested_xml)
    else: # parallel_process
        for chunk in chunks(query_generator, n_compute_nodes):
            # mapper() applies process_query to each element of chunk
            # this is where parallel processing occurs.
            # "chunk" is a list of n_compute_nodes queries
            try:
                for (xmlok, query, requested_xml), output in pool.map(process_query, chunk):
                    log.info("BIT Query:", repr(query))
                    log.write(output, eol="")
                    if xmlok:
                        simple_file.write(requested_xml)
                    if not xmlok or "<message>" in requested_xml:
                        log.error("Exception:", requested_xml)
            except KeyboardInterrupt as e:
                pool.terminate()
                raise

    simple_file.write("</Responses>\n")

def chunks(generator, chunk_len):
    """Chop generator of large sequence into smaller list chunks since
    IPython parallel processing doesn't work on generators, just
    things with __len__().  Also yield multiple smaller sequences to
    avoid pre-computing 70k queries.
    """
    chunk = []
    for item in generator:
        chunk.append(item)
        if len(chunk) >= chunk_len:
            yield chunk
            chunk = []
    yield chunk

# ========================================================================

def query_to_peng(query_str):
    """ Translate a single query string to pyetc engine inputs (peng) dict. """
    web_inputs, requests = translate.process_query(query_str, info=False)
    engine_inputs = compute.web_to_engine_inputs(web_inputs)
    return web_inputs, engine_inputs

def debug_queries(query_list):
    """ Loop over them. """
    for q in query_list:
        print('\n'+q)
        debug_query(q)
        junk = input('Okay?')

def run_query(query_str, detectors, only_response=False, debug=False):
    """Process and dump a single query or a file of queries.

    Parameters
    ----------
    query_str: string
        The query to process.
    detectors: list
        The list of HST detectors ('ccd', 'nuvmama', etc.)
    only_response: bool
        If False, this function will return a complete XML BIT table (including versions) for the query.
        If True, this function will only return the <Response>...</Response> element for the query.

    Returns
    -------
    xml_response: string
        The XML response for the query, formatted as specified by the only_response parameter.
    """
    log.VERBOSE_FLAG = True

    web_inputs, engine_inputs = query_to_peng(query_str)

    # HERE we call peng3jeng
    jeng_dict = peng3jeng.pyetc_in_dict_to_pandeia_in_dict(engine_inputs, None, pweb_dict=web_inputs)

    if debug:
        log.verbose("\npandeia engine inputs:\n", pprint.pformat(jeng_dict), "\n")

    # Run the Pandeia calculation on the jeng input.
    jeng_results = perform_calculation(jeng_dict, webapp=True, dict_report=True)

    if debug:
        log.verbose("\npandeia engine outputs:\n", pprint.pformat(jeng_results), "\n")

    # Get the versions string.
    versions_str = get_version_string()

    # Extract what we need from the results and convert to XML.
    xml_response = ''
    if not only_response:
        xml_response = '<Responses>\n'
        xml_response += versions_str + '\n'

    xml_response += create_xml_response(query_str, jeng_results, detectors) + '\n'

    if not only_response:
        xml_response += '</Responses>'

    if debug:
        log.info("XML:", xml_response)

    return xml_response

# ========================================================================


def main():
    """Process a file of bit queries, `input`,  into a file of XML results,
    `output`.
    """
    parser = optparse.OptionParser(usage="usage: %prog [options] <input> <output>")

#   parser.add_option(
#       "-n", "--parallel-process-nodes", dest="parallel_compute_nodes",
#       help="Specify CPU count for parallel processing,  default=1.",
#       metavar="NODES", default=1)

    parser.add_option(
        "-q", "--query", dest="query",
        help="URL query string for a BIT request",
        metavar="QUERY", default=None)

    parser.add_option("-f", "--queryfile", dest="query_file", help="Input file with one query per line", default=None)

    parser.add_option("-o", "--outdirfn", dest="output_dir_filename", help="Full path and filename of output table", default="./bit_table.xml")

    parser.add_option('-D', '--direct', dest='direct', action='store_true',
                      help="Direct call - do not employ pdb wrapper")

#   options, args = log.handle_standard_options(parser, sys.argv, show_info=False)
    options, args = parser.parse_args(sys.argv)

    # Get the list of detectors.
    detectors = get_detectors_list()

    if options.query_file is not None:
        # A query file is set, so read the queries from there.
        query_file = open(options.query_file, 'r')
        queries = query_file.readlines()
        query_file.close()

        xml_responses = '<Responses>\n'
        xml_responses += get_version_string()

        # Process each of the queries.
        for query in queries:
            try:
                response = run_query(query, detectors, only_response=True)
                xml_responses += response
            except Exception as e:
                log.info(e)

        xml_responses += '</Responses>'

        out_file = open(options.output_dir_filename, 'w')
        out_file.write(xml_responses)
        out_file.close()
    elif 0 and options.query is None:
        # cds: TRANSLATE THIS to get rid of standard_run()
        log.standard_run("process_queries(args[1], args[2], "
                         "int(options.parallel_compute_nodes))",
                         options, globals(), locals())
    else:
        if 1 or options.direct: # run ONLY direct for a hot minute to get this started
            try:
                response = run_query(options.query, detectors)
                out_file = open(options.output_dir_filename, 'w')
                out_file.write(response)
                out_file.close()
            except Exception as e:
                log.info(e)
        else:
            log.standard_run("run_query(options.query)",
                             options, globals(), locals())
    log.write()
#   log.standard_status()

#
# For example, to run a single query (to debug say), it is
# % ./bit.py -D -q 'type=bot&amp;magnitude=10&amp;config=.....'
#
# Or, to generate a table from a file containing multiple queries, use:
#  python bit.py -f /path/to/query_file.txt -o /path/to/output_directory/output_table.xml
#  bit_table.xml will be written out to the specified directory.
if __name__ == "__main__":
    main()
