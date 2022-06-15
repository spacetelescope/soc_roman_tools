'''Map the extraction form parameters into a form that the etc_engine can use.
'''

from __future__ import division

import os

#from ...main import util, tools
from pandeia.engine.helpers.bit.pyetc_util import read_dict
import pandeia.engine.helpers.bit.instruments as instruments


mapfile = os.path.join(os.path.dirname(__file__), 'etc_web_data/extraction.dat')
EXTRACTION_MAP = read_dict(mapfile)


def _handle_scan_mode(ans, fd, shape):

    # Here, we handle the 'width' value.
    #
    # 'varname' is the name of the HTML field that contains the relevant extraction parameter:
    # - width for imaging. The height is computed elsewhere.
    #   - height (scan length) for spectroscopy. The width is implicitly assumed to be 1 resel.
    #
    # The logic here assumes that the scan direction is parallel to the "height" direction.
    # We have to handle the fact that the extraction rectangle dimensions in spectroscopy are
    # reversed wrt imaging (an engine API hurdle).
    if fd['science_mode'] == 'scimaging':
        # in imaging, the dimension that matters for PSF-extraction purposes is the width (perpendicular
        # to the scan direction). It is used to index and interpolate the EE table. It is directly read
        # from the page. It will end up generating a rectangle spec such as [5.0, 82.57]. In here, the
        # first element (5.0) is the width read from the page, and the second element (82.57) is computed
        # elsewhere.
        varname = "extractionRegionScWidth"

    else:
        # in spectroscopy, there is no dimension associated with EE table weighting. The width is implicitly
        # 1 resel, and it is inserted into the extraction region specification later on. The only relevant size
        # is the scan length itself, which is the height of the extraction region. The same scan length used
        # in the above comment would then generate a rectangle spec like [82.57, 1.0], meaning 82.57 pixels
        # high and 1 resel wide.
        varname = "extractionRegionScHeight"

    # Look up the size (in pixels) from the variable name.
    size = fd[varname]
    ans['size'] = size

    # Scan mode adds further complexity to rectangles. Here, we handle the 'height' value.
    if 'rectangle' in shape:

        height = fd['extractionRegionScHeight']

        # 'extractionRegionScHeight' is a drop-down selector
        # that has entries of 'scanLength' (assumed to be in
        # arcsec) and '1' (assumed to be in pixels).

        # Only run the clause in case the drop-down was selected
        # as 'scanLength'. The '1' case is handled trivially by
        # just holding the value, unaltered, in the 'height' variable.
        if height == 'scanLength':

            # scan length can be gotten from either of two page variables
            # Which one to use is controlled by the 'simmode' radio button.
            # This radio button can have the values of either 'SNR', or
            # 'SNR_from_rate'. We *have* to use 'SNR' for the button case
            # that allows the user to enter the Time value, so as to re-use
            # the existing web code logic for stare mode.

            # 'SNR_from_rate' means: "compute SNR given a Time value derived
            # from a scan rate and a scan length".
            if fd['simmode'] == 'SNR_from_rate':
                height = fd['scLengthFromRate']
            # 'SNR' means what it always meant in stare mode: "compute SNR
            # given a Time derived from a page field"
            if fd['simmode'] == 'SNR':
                height = fd['scLengthFromTime']

        # need plate scale to convert extraction height to pixels.
        instrument_name = fd['instrument']
        instconfig = tools.find_instrument('hst', instrument_name)
        dconfig = instconfig.data.get_dconfig(fd['detector'])
        pixel_height = dconfig['pixel_height']

        # in scan mode, when height=1, it really means height=1 pixel.
        # Any other value, it's assumed height is in arcsec and thus
        # must be converted to pixels before feeding to the engine.
        height = float(height)
        if height != 1.0:
            height /= pixel_height

        # Extraction box in imaging is reversed wrt spectroscopy...
        if fd['science_mode'] == 'scimaging':
            ans['size'] = [float(size), float(height)]
        else:
            ans['size'] = [float(height), 1.]


#TODO web cleanup:
# this function is too long; should be broken into smaller pieces. Obvious candidate is the scan-mode-specific code.

def extraction_info(formdata, scan_mode):
    """ Map form data into extraction data, with instrument-specific
    lookups to provide 'default' (aka recommended) values.
    """
    fd = formdata
    ans = {}

    # Copy in the gyromode (for now at least)
    ans['gyromode'] = fd['gyromode']

    # Look up the source type (point or extended)
    target_type = fd['fsourceType']

    # Based on that, look up the source diameter
    # (in arcsec) and the aperture shape.
    shape = _handle_target_type(ans, fd, target_type)

    #TODO web cleanup:
    # probably need to add a hook for instrument-specific customizations here

    # BUG: in the meantime there is a hardcoded hack for
    # COS specacq mode which needs an extra parameter
    if (fd['instrument'] == 'cos' and
        fd['science_mode'] == 'spectroscopicacq'):
        # then we need an extra variable
        # (only for NUV, because it makes no difference for FUV)
        ans['acq_mode'] = fd.get('nuvMode','')
        ans['stripe'] = fd.get('Stripe','')

    # (Source type, aperture shape) key the extraction_map, which
    # determines which variable to look in. Disperser is used in
    # a subset of cases.
    disperser_ = fd.get('disperser', '')

    if shape == 'default':
        # Then look up the instrument-specific defaults, which
        # depend on source type and mode as well as detector.
        ans['shape'],ans['size'] = _get_inst_extraction_defaults(fd['instrument'],
                                                                 fd['science_mode'],
                                                                 fd['detector'],
                                                                 disperser_,
                                                                 target_type)
    else:
        # Otherwise, read it from the standard map
        ans['shape'], varname = EXTRACTION_MAP[(target_type, shape)]

        # there might be one level of indirection. This happens for
        # now only with wfc3 uvis and ir spectroscopic mode.
        if varname == 'xRegionType':
            rt = fd[varname]
            if rt == 'default':
                ans['shape'],ans['size'] = _get_inst_extraction_defaults(fd['instrument'],
                                                                         fd['science_mode'],
                                                                         fd['detector'],
                                                                         disperser_,
                                                                         target_type)
                # invalidate so it won't be used after here.
                varname = ''
            else:
                ans['shape'], varname = EXTRACTION_MAP[(target_type, rt.lower())]

        # Look up the size (in pixels) from the variable name. This will work for
        # most variable names, but fails with the special variable names used in
        # scan mode. Scan mode needs special handling.
        if varname in fd:
            size = fd[varname]
            ans['size'] = size

            # spectroscopic pages deliver the extraction aperture size as a
            # comma-separated string. Convert it to the engine API format, a
            # tuple. Conversions look unsafe but will be done only on values
            # provided by the html code itself. In case one envisions other
            # uses for this code, a safer version must be put in place.
            if type(size) == type(""):
                size = size.split(",")
                ans['size'] = (float(size[0]), float(size[1]))

        # Scan mode adds complexity to input page fields.
        if scan_mode:
            _handle_scan_mode(ans, fd, shape)

    #TODO web cleanup:
    # an instrument hook here would be handy.

    if (fd['instrument'] == 'stis' and
        fd['science_mode'] == 'targetacquisition' and
        "peak" in fd['ccdMode'].lower()):

        ans['shape'] = 'square'
        ans['size']  = 32

    # #815: IMAGE / SEARCH modes for COS image target acquisition.
    if (fd['instrument'] == 'cos' and
        fd['science_mode'] == 'targetacquisition'):

        if fd.get('cosAcqMode','') == 'search':
            ans['shape'] = 'rectangle'
            ans['size']  = (220,470)

    #That's it.
    return ans


def _handle_target_type(ans, fd, target_type):

    ans['stype'] = target_type

    if target_type == 'point':
        ans['sdiameter'] = 0.
        shape = fd.get('xRegionType', '').lower()

    elif target_type == 'extended':
        ans['sdiameter'] = fd['fdiameter']  # fdiameter is in arcsec
        shape = fd.get('xRegionExtendedType', '').lower()

    else:
        raise ValueError('Unknown type of extraction region')

    return shape


def _get_inst_extraction_defaults(instrument,
                                  mode,
                                  detector,
                                  disperser,
                                  target_type): #point/extended



    """ Return info from instrument-specific dictionary.
    Be smart and cache results so we don't have to repeatedly read the same
    dictionary file.

    """
    # TODO:
    #It's possible this should call out to a module in the specific
    #instrument, because the default dictionaries are still being shoehorned
    #into a standard.

    # Import the correct module
    mod = instruments.import_instrument_module(instrument,
                                                   '_web')
    # (Efficiently) get the data we need
    EXTRACTION_DEFAULTS = mod.get_data('extraction_defaults')

    # Now it has the right dictionary loaded: look up the
    # default values.
    # Some instruments may have the default extraction regions
    # keyed with 3 parameters, most are keyed with 2.
    try:
        defval = EXTRACTION_DEFAULTS[(target_type, mode)]
    except KeyError:
        index = disperser.find('_0')
        grating = disperser[index+2:]
        defval = EXTRACTION_DEFAULTS[(target_type, mode, grating)]

    # Some instruments have a further level of nesting by dict;
    # others are simpler.
    if isinstance(defval, dict):
        shape,size = defval[detector]
    else:
        shape, size = defval

    return shape, size
