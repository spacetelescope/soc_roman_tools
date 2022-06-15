import os

import pandeia.engine.helpers.bit.pyetc_util as util
util.INST_DATADIR = os.path.dirname(__file__)

from pandeia.engine.helpers.bit.instruments.telescopes import TELESCOPES
#####

def import_instrument_module( instrument, module='' ) :
    '''
    import 'instruments.hst.%s.%s%s'%(name,name,which)

    for example:
       engine_module = import_instrument_module( 'acs', '' )
       web_module    = import_instrument_module( 'acs', '_web' )

    '''
    # instrument is the name of the instrument that we are interested in.
    # After #1170, web code lives under a 'web' subpackage. Note also that
    # after the engine refactoring, this method is only called with a
    # module='_web' option. Unit tests call it with module='', and nobody
    # else uses the option.
    if module == '_web':
        name = 'pandeia.engine.helpers.bit.instruments.hst.%s.web.%s%s' % ( instrument, instrument, module )
    else:
        name = 'pandeia.engine.helpers.bit.instruments.hst.%s%s' % ( instrument, module )

#   try :
    return util.dynamic_import(name) # for BIT runs we can handle if this crashes here to highlight any problem
#   except ImportError:
#       pass

    # We could try other possible names here in a similar manner.  When users
    # can provide their own instrument package, it will be convenient to
    # allow them to put the package someplace other than pyetc.instruments.hst
    #

    # if we get to the end, re-raise the last import error
    raise


#####

# This is a short-lived cache for get_instrument_data; it allows us to
# get_instrument_data() whenever we want the data, but we don't have to
# actually read the file multiple times.
instrument_data_cache = { }

def get_instrument_data( instrument, name ) :
    '''
    Read an instrument-specific data file.
    (replaces most uses of etc_engine.util.read_dict())

    For example:
        d = instruments.get_instrument_data( 'acs', 'sharpness' )

    This capability is also present under the name get_data() in each
    instrument-specific web module:

        m = instruments.import_instrument_module( 'acs', '_web' )
        d = m.get_data( 'sharpness' )

    In either case, the file is only read once; the same copy of the
    file data is returned each time.  (Take care not to alter the
    returned object.)

    '''
    index = ( instrument, name )

    if index in instrument_data_cache :
        return instrument_data_cache[index]

    # Find the data directory for this instrument
    dirname = os.path.join(util.INST_DATADIR, 'hst', instrument)

    # read the data
    fname = '%s/%s_%s.dat' % ( dirname, instrument, name )

    datadict = util.read_dict(fname)

    # cache it and return it
    instrument_data_cache[index] = datadict

    return datadict

def get_alt_instrument_data( instrument, name ) :
    '''
    Read an instrument-specific data file.
    (replaces most uses of etc_engine.util.read_dict())

    For example:
        d = instruments.get_instrument_data( 'acs', 'sharpness' )

    This capability is also present under the name get_data() in each
    instrument-specific web module:

        m = instruments.import_instrument_module( 'acs', '_web' )
        d = m.get_data( 'sharpness' )

    In either case, the file is only read once; the same copy of the
    file data is returned each time.  (Take care not to alter the
    returned object.)

    This function should be called when there is a possibility that
    an alternate version of the same data file exists in the
    telescope/instruments/ sub-directory. For now we must code this
    as two separate functions since this module is used by both the
    old and the new versions of the engine code.

    '''
    index = ( instrument, name )

    if index in instrument_data_cache :
        return instrument_data_cache[index]

    # Find the data directory for this instrument
    dirname = os.path.join(util.INST_DATADIR, instrument)

    # read the data
    fname = '%s/%s_%s.dat' % ( dirname, instrument, name )

    datadict = util.read_dict(fname)

    # There might be a new version of he file in pyetc/instruments/.
    # Read it and override former version just read.
    new_name = fname.replace('instruments/hst/%s'%instrument,
                             'instruments/hst/%s/%s/'%(TELESCOPES[instrument],instrument))
    try:
        datadict = util.read_dict(new_name)
    except IOError as e:
        pass

    # cache it and return it
    instrument_data_cache[index] = datadict

    return datadict

#####

# Wrapper to the specific
# <instr>_obsmode.make_<instr>_<sciencemode>_obsmode function

def make_obsmode_test_interface(form_input):
    """This is a test rig to find a back door into instrument specific
    obsmode functions.  For non-test code, you would call the correct
    function directly from the relevant instrument-specific code.
    """

    instr = form_input['instrument'].lower()
    mode = form_input['science_mode'].lower()

    m = import_instrument_module(instr, '_web')

    funcname = 'make_' + instr + '_' + mode + '_obsmode'
    if funcname in m.__dict__ :
        fn = m.__dict__[funcname]
        ans = fn(form_input)
        return ans

    raise NotImplementedError('Instrument/mode combination %s, %s not supported'%(instr,mode))

#####

#Wrapper to the specific <instr>.instrument_from_form function
def instrument_from_form(instr_name, science_mode, form_input):
    """Call the instrument-specific function to process the web form inputs
    into the instrument-specific dictionary.

    instr_name:   lower case string
    science_mode: lower case string
    form_input:   dictionary containing form inputs in proper type
    """

    #This is specific to each instrument_sciencemode combination
    instr_namel = instr_name.lower()

    m = import_instrument_module(instr_namel, '_web')
    idict = m.instrument_from_form(science_mode, form_input)
    return idict, idict.pop("wodict", {})

#####

# This list of all instruments is here for testing only.
test_mode_all_instruments = (
    'acs',
    'cos',
    'stis',
    'wfc3ir',
    'wfc3uvis',
)
