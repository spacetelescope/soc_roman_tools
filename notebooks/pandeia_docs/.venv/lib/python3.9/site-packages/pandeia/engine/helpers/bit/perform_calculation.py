from __future__ import division
import errno, re, sys

from pandeia.engine.helpers.bit.pyetc_util import read_dict, PYSYN_CDBS
from pandeia.engine.helpers.bit import instruments
from pandeia.engine.helpers.bit import string_constants as sc
from pandeia.engine.helpers.bit.extraction import extraction_info
from pandeia.engine.helpers.bit.sky import Sky
from pandeia.engine.helpers.bit.source import Source

# hard-coding subsystem since this code is ONLY here in pandeia to run the BIT
#subsystem = 'etc'
subsystem = 'bit'

# We need to know what telescope to use to find the MJD for the current
# calculation.  This version of the ETC does not really understand multiple
# telescopes, so we assume that we are only configured with one.
telescope = 'hst'

# This is the MJD dict, lazily loaded from the $CDBS/etc/mjd.dat
mjd_dict = None

# obtain the MJD to use for a particular configuration
#
# subsystem normally comes from etc_web.config.subsystem
#
# telescope is which telescope to use.  Since etc_web doesn't currently
# understand multiple telescopes, use the global value telescope from
# this file.
#
# instrument is the name of the instrument.
#
def _handle_MJD( subsystem, telescope, instrument, engine_input ):

    # lazy load mjd_dict
    global mjd_dict
    if mjd_dict is None :
        try :
            mjd_dict = read_dict( PYSYN_CDBS + '/etc/mjd.dat' )
        except IOError as e :
            if e.errno != errno.ENOENT :
                raise
            # no file means there are no entries
            mjd_dict = { }

    # The priorities are:
    # 1. a fully specified entry
    # 2. a subsystem/telescope entry
    # 3. None
    try :
        mjd = mjd_dict[ ( subsystem, telescope, instrument ) ]
    except KeyError :
        try :
            mjd = mjd_dict[ ( subsystem, telescope ) ]
        except KeyError :
            return

    # try :
    #     mjd = mjd_dict[ ( subsystem, telescope, instrument ) ]
    # except KeyError, e :
    #     raise exceptions.DataError("Invalid entry for MJD calculation. Exception is: " + str(e))

    # add MJD to obsmode
    engine_input['obsmode'] += ',MJD#' + str(mjd)


#TODO web cleanup:
#  this function is too long; does too many things.

def process_form_input(instrument, science_mode, form_input, request_dir):
    """ Process the form input dictionary to produce the inputs that the
    engine wants. Return the engine input."""

    # web_output is computation performed here in the web code
    web_output = {}

    # engine input is input to be passed on to the engine
    engine_input = dict()

    ## bug: Apparently, these should be in all of the input forms, but are not.
    #TODO web cleanup:
    #  so, fix it! make sure all forms have 'gyromode', and remove this statement.
    if not 'gyromode' in form_input:
        form_input['gyromode'] = 'normal'

    #and the sky expression & descriptions
    sky = Sky(formdata=form_input)
    sky.make_synexpr()
    engine_input['sky_expr'] = sky.expr
    engine_input.update(sky.descr)

    # instrument_from_form() calls into the instrument
    # specific instrument_from_form(), which converts selected fields
    # from the form input into a dictionary of values
    # to pass on to the engine.
    #
    # obsmode is one of those values.
    idict, wodict = instruments.instrument_from_form(instrument,
                                                     science_mode,
                                                     form_input)
    engine_input.update(idict)
    web_output.update(wodict)

    # add MJD to the instrument obsmode
    _handle_MJD( subsystem, telescope, instrument, engine_input )

    # Create the source expression. This must be created after the
    # instrument-related values have been populated in order to support
    # normalization in the observation bandpass.

    # Due to a bug in the pysynphot parser, the obsmode has to be
    # validated against problematic aperture specifications. Do that first.
    normalization_obsmode = _validate_obsmode(engine_input.get('aperture'),
                                              engine_input.get('obsmode'),
                                              form_input.get('fftype')
                                              )

    src = Source(request_dir=request_dir, formdata=form_input,
                 obsmode=normalization_obsmode)

    # Use this to populate the relevant engine_input and web_output items.
    engine_input['source_expr'] = src.expr
    engine_input['target_lines'] = src.lines
    web_output['source_help'] = src.help

    #Create the extraction info

    scan_mode = science_mode in [sc.SCIMAGING, sc.SCSPECTROSCOPIC]

    # The extraction module looks for the detector in form_input
    #TODO web cleanup:
    # put it in all the forms, then remove this hack
    form_input['detector'] = engine_input['detector']
    engine_input['extraction'] = extraction_info(formdata=form_input, scan_mode=scan_mode)

    #The extraction info wants the detector: add it.
    #TODO web cleanup:
    # couldn't this be done by extraction_info itself ?
    try:
        engine_input['extraction']['detector'] = engine_input['detector']
    except KeyError:
        pass

    #Treat some specific elements in special ways
    engine_input['instrument'] = instrument
    engine_input['science_mode'] = science_mode

    if "Time" in form_input:
        # Normal operation.
        engine_input['time'] = form_input['Time']

    if not scan_mode:
        #If the keys are spelled the same, we can just copy them
        klist = ['SNR',
                 'simmode']
        for key in klist:
            engine_input[key] = form_input[key]

        if 'spectroscopic' in science_mode:
            web_output['spectroscopic'] = True
        else:
            web_output['imaging'] = True

        web_output['mode'] = science_mode

    else:
        _handle_scan_mode(engine_input, form_input, web_output)

    # Rampfilter needs to be set to imaging;
    # otherwise use the science_mode as is.
    if science_mode == 'rampfilter':
        web_output['mode'] = 'imaging'

    return engine_input, web_output


def _handle_scan_mode(engine_input, form_input, web_output):
    """ Handle the specifics associated with scan mode.

    Scan mode supports two possible user inputs:
      - user enters scan rate and scan length. The code in here
        computes Time, which is then fed to the ETC engine to
        compute SNR.
      - user enters Time and scan length. Time is fed to the ETC
        engine as above to compute SNR. The code in here computes
        the scan rate for reporting and warning purposes only.

    The input dictionaries are modified in place.

    Parameters
    ----------
    engine_input: dict
      dict with the engine input
    form_input: dict
      dict with the form input associated with the web request
    web_output: dict
      dict with values used by the output templates

    """
    # Get scan rate from the input request. This might be
    # overwritten below, depending on other input settings.
    engine_input['scrate'] = form_input['scRate']

    # There are two possible sources of scan length in the input
    # page. Pick the one selected by the 'simmode' radio button.
    key = 'scLengthFromRate'
    if form_input['simmode'] == 'SNR':
        key = 'scLengthFromTime'
    engine_input['sclength'] = form_input[key]

    # Time is computed once and for all from the scan length and rate.
    # Or, it's just read from the Time field.
    # Page and/or validation code should handle limiting cases.
    simmode = form_input['simmode']
    if simmode == 'SNR_from_rate':
        engine_input['time'] = float(engine_input['sclength']) / float(engine_input['scrate'])
    else:
        engine_input['time'] = form_input['Time']
        engine_input['scrate'] = engine_input['sclength'] / engine_input['time']

    # Time is reported directly. No inverse computation is supported.
    web_output['time'] = engine_input['time']

    # SNR is passed to the engine input but never used in scan mode.
    # This only fulfills the engine API requirements.
    engine_input['SNR'] = form_input['SNR']

    # We only support the direct SNR calculation.
    engine_input['simmode'] = 'SNR'

    # report correct mode at output.
    if 'scspectroscopic' in engine_input['science_mode']:
        web_output['scspectroscopic'] = True
    else:
        web_output['scimaging'] = True

    # Some science modes look better when reported
    # in a more human-readable form.
    if engine_input['science_mode'] == 'scimaging':
        web_output['mode'] = 'imaging scan'
    elif engine_input['science_mode'] == 'scspectroscopic':
        web_output['mode'] = 'spectroscopic scan'


def _validate_obsmode(aperture, obsmode, normtype):
    """ Validate for avoidance of pysynphot parser bug (pysyn #238)

    Some apertures/slits are spelled with digits followed by letters,
    and this causes the pysynphot parser to incorrectly break up the
    keyword which then causes an exception either in the parsing phase
    or in the graph table lookup. This function removes problematic
    apertures from the obsmode string.

    Parameters
    ----------
    aperture: string
    obsmode: string
    normtype: string

    Returns
    -------
    result: string

    """
    result = obsmode
    # The bug appears only when normalizing by bandpass.
    if normtype != 'fNormalizeByBandpass':
        return result

    # Not all obsmodes have apertures.
    if aperture is None:
        return result
    if aperture not in obsmode:
        return result

    # Not all apertures are problems.
    # Use re.match to determine whether it starts with a digit.
    problem = re.match("[0-9]", aperture)
    if problem:
        raise ValueError("Normalization in observation bandpass that includes aperture %s is not supported at this time."%aperture)

    return result
