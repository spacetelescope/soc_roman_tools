import pandeia.engine.helpers.bit.instruments

#
# read a dictionary from one of the .dat files

def get_data( name ) :
    return pandeia.engine.helpers.bit.instruments.get_instrument_data( 'acs', name )

FILTERS = get_data("config")["filters"]
DETECTOR_TYPE = dict(hrc='ccd',wfc='ccd',sbc='mama')

def instrument_from_form(science_mode, form_input):
    """Function to extract values of interest from the web form
    data and return a smaller dictionary which will be used to
    update the engine_input dictionary.
    """
    idict = dict( gain =   form_input['gain'],
                  nreads = form_input['crsplit'],
                  fdiameter = form_input['fdiameter'],
                  )
    idict["wodict"] = wodict = {}
    wodict["coronography"] = "off"

    #check if nreads is zero and return an error
    if idict["nreads"]  < 1:
        raise ValueError("CRSPLIT cannot be less than 1")

    if science_mode == 'imaging':
        #make an obsmode string that can be submitted to pysynphot
        idict['obsmode'], filter, polarizer = make_acs_imaging_obsmode(form_input)
        idict['detector'] = form_input['detector']
        detector = idict['detector']
        if idict['detector'].lower() == 'sbc':
            idict['gain']=1.
        idict["filter"] = wodict["filter"] = filter
        if filter is not None:
            idict["filter"] = filter.lower()
        elif polarizer is not None:
            # this should be a case of clear filter plus polarizer so use
            # polarizer for bandpass limits.
            idict["filter"] = polarizer.lower()
        else:
            # this should be the case where there are two clear filters so pick one.
            if detector == 'wfc':
                idict['filter']=form_input['wfcfilt1'].lower()
            elif detector == 'hrc':
                idict['filter']=form_input['hrcfilt1'].lower()
        wodict["filter"] = idict["filter"]
        if filter is not None:
            wodict["filter_info"] = FILTERS[filter]
    elif science_mode == 'rampfilter':
        #make an obsmode string that can be submitted to pysynphot
        idict['obsmode']=make_acs_rampfilter_obsmode(form_input)
        wodict["rampfilter"] = True
        detector=form_input['detector']
        idict['detector'] = detector
        if detector == 'wfc':
            idict['filter']=form_input['wfcfilt1'].lower()
        elif detector == 'hrc':
            idict['filter']=form_input['hrcfilt1'].lower()
        wodict["filter"] = idict["filter"]
        wodict["filter_info"] = FILTERS[idict["filter"]]
        polarizer = None
    elif science_mode == 'spectroscopic':
        #the spectroscopic form is set up so that the user chooses a disperser
        #and the name of the disperser contains the name of the detector. I don't
        #like this, but am keeping it for now.
        combine_string=(form_input['disperser']).split('_')
        detector=combine_string[0].lower()
        grating=combine_string[1].lower()[1:]
        idict['obswave']=form_input['obswave']
        idict['grating']=grating
        idict['detector']=detector

        if detector == 'wfc':
            filter=form_input['wfcfilt1'].lower()
        elif detector == 'hrc':
            selector = {"0":"1", "1":"0"}[combine_string[1][0]]
            filter=form_input['hrcfilt' + selector].lower()
            wodict["filter"] = filter
        elif detector == 'sbc':
            filter='clear' #to pass the concatination
            idict['gain']=1
        idict["filter"] = filter
        #make an obsmode string that can be submitted to pysynphot
        idict['obsmode'] = make_acs_spectroscopic_obsmode(idict)

        wodict["disperser_info"] = disperser = {}
        disperser["label"] = {"p" : "Prism", "g" : "Grism"}[grating[0]]
        disperser["grating"] = grating.upper()
        disperser["name"] = FILTERS[grating]["name"]
        if not filter.startswith("clear"):
            polarizer = filter
        else:
            polarizer = None
    else:
        raise NotImplementedError('%s mode not implemented for acs' % science_mode)
    if polarizer:
        wodict["polarizer"] = {
            "pol_v":"Visible",
            "pol_uv":"Ultraviolet"
        }[polarizer]

    wodict['detector_type'] = DETECTOR_TYPE[idict['detector']]

    # Silently discard post_flash inputs for non-CCDs and non-supported spectrograph mode
    # See #922
    wodict["show_post_flash"] = False
    # NOTE that w/ show_post_flash, we add a control-flow flag in a place
    # where its not expected. During code review jwcalibdev.stsci.edu/r/69
    # we decided to just comment for now.  If it happens again, we may want
    # to add a dict under/near wodict to handle these things.
    if science_mode != 'spectroscopic' and \
       DETECTOR_TYPE[idict['detector']] == 'ccd':
        idict['post_flash_electrons'] = form_input['post_flash_acs']
        wodict["show_post_flash"] = True

    return idict

def make_acs_imaging_obsmode(form_input):
    # instrument
    elements = ['acs']
    #The detector gets an entry of its own on only some of the forms
    #on other forms it derived from the value of one of the other
    #varibles (as in spec  mode). If that's the case pass in the detector
    #to this function.
    detector = form_input['detector']
    #QUERY: wfc1 and wfc2 are distinct obsmodes,
    #representing the two halves of the chip. But the
    #ETC always assumes wfc1.
    if detector == 'wfc':
        elements.append('wfc1')
    else:
        elements.append(detector)
    #filters: Based on the detector, we decide which filters to pick out.
    f0 = form_input['%sfilt0' % detector].lower()
    if detector != 'sbc':
        #there are two filter wheels, either of which may be clear
        f1 = form_input['%sfilt1' % detector].lower()
        filters = [f for f in (f0, f1) if not f.startswith('clear')]
        elements.extend(filters)
    else:
        #there is only one and it can't be clear
        elements.append(f0)
        filters = [f0]
    obsmode = ','.join(elements)

    # For now anyway,  polarizers appear on the form with a phase but are
    # only uploaded to the server as V or UV.
    polarizer = filter = None
    filter_count = 0
    for f in filters:
        if f.startswith("pol"):
            polarizer = f
        elif not f.startswith("clear"):
            filter = f
            filter_count += 1

    #Enforce restriction against crossed filters
    if filter_count > 1:
        raise ValueError("Use of two non-clear filters %s is not supported. Please select only one filter and an optional polarizer." % str(filters))

    return obsmode, filter, polarizer

def make_acs_rampfilter_obsmode(form_input):
    #we need a string like: acs,wfc,fr388n#3880

    #instrument
    elements = ['acs']
    #Detector
    detector = form_input['detector']
    #wfc always uses wfc1
    if detector == 'wfc':
        elements.append('wfc1')
    else:
        elements.append(detector)

    #filter: choose based on detector
    f = form_input['%sfilt1'%detector].lower()
    #and the wavelength for the ramp filter
    wave=form_input['obswave']

    elements.append("%s#%s"%(f,wave))
    obsmode = ','.join(elements)

    # There is a KNOWN BUG in pysynphot right now that
    # does not support the full range of wavelengths for the
    # ramp filters. So try to construct a bandpass now, so
    # we can raise an appropriate exception if it fails.
    #
    # This is problematic because it requires importing
    # pysynphot into this web module, which in general is
    # disallowed for clean-interface reasons. This code
    # should be removed as soon as the known pysynphot bug
    # is fixed.
    #
    # TODO: remove this code when pysynphot is fixed
    # TODO: think about better error message management.
    import pysynphot as psyn
    try:
        psyn.ObsBandpass(obsmode)
    except psyn.exceptions.ExtrapolationNotAllowed as e:
        msg = """You have selected instrument configuration %s at a wavelength of %g. This wavelength is outside the range
              of this filter. Please select a wavelength that is within the supported range for this filter."""
        raise ValueError(msg%(str(elements[:]), wave))
    return obsmode

def make_acs_spectroscopic_obsmode(form_input):
    #instrument
    elements = ['acs']

    # This function can be called with a 'processed' input dictionary,
    # as well as a raw one. The raw dictionary contains exactly what
    # is in the input page. The input spectroscopic page has no 'detector'
    # 'grating', or 'filter' fields, so we get them from the 'disperser'
    # string. This is similar to what is done in instrument_from_form()
    # above, but not quite the same. An effort should be put here when
    # refactoring the web code.
    if 'detector' in form_input:
        detector = form_input['detector']
        grating = form_input['grating']
        filter = form_input['filter']
    else:
        combine_string = (form_input['disperser']).split('_')
        detector = combine_string[0].lower()
        grating = combine_string[1].lower()[1:]
        if detector == 'wfc':
            filter = form_input['wfcfilt1'].lower()
        elif detector == 'hrc':
            selector = {"0":"1", "1":"0"}[combine_string[1][0]]
            filter = form_input['hrcfilt' + selector].lower()
        elif detector == 'sbc':
            filter = 'clear'
        else:
            filter = ''

    #detector: wfc always uses wfc1
    if detector == 'wfc':
        elements.append('wfc1')
    else:
        elements.append(detector)

    #Based on the detector, we decide which filters to pick out.
    f = [x for x in (grating, filter) if not x.startswith("clear")]
    elements.extend(f)
    obsmode = ','.join(elements)
    return obsmode
