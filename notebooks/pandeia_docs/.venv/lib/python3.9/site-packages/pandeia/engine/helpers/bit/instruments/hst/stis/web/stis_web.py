from __future__ import division

import pandeia.engine.helpers.bit.instruments

DETECTOR_TYPES = {'ccd' : 'ccd', 'fuvmama': 'mama', 'nuvmama': 'mama'}

#
# read a dictionary from one of the .dat files
def get_data( name ) :
    return pandeia.engine.helpers.bit.instruments.get_instrument_data( 'stis', name )

GRISMS = get_data("config")["gratings"]
CI_TABLE = get_data("config")["ci_table"]


def instrument_from_form(science_mode, form_input):
    """Function to extract values of interest from the web form
    data and return a smaller dictionary which will be used to
    update the engine_input dictionary.
    """
    # This gets rid of certain complaints from the migration code.
    binx    = form_input.get("binx",1) #get_default(form_input, "binx")
    biny    = form_input.get("biny",1)
    crsplit = form_input.get("crsplit",1)
    gain    = form_input.get("gain",1)
    #TODO: refactor so some common values can be set upfront
    if science_mode == 'spectroscopic':
        obsmode = make_stis_spectroscopic_obsmode(form_input)
        disperser, detector, grating, central_wavelength, aperture = _get_stis_spec_parameters(form_input)
        idict = dict(disperser=disperser,
                     detector=detector,
                     grating=grating,
                     aperture=aperture,
                     obsmode=obsmode,
                     crsplit=crsplit,
                     gain=gain,
                     binx=binx,
                     biny=biny,
                     )
        idict["wodict"] = wodict = {}

        idict['obswave'] = form_input['obswave']
        if central_wavelength != None:
            idict['central_wavelength'] = central_wavelength

        wodict["grating_info"] = GRISMS[grating.lower()]
        wodict["aperture"] = obsmode.split(',')[-1]
        wodict["detector_type"] =  DETECTOR_TYPES[detector].lower()

    elif science_mode == 'imaging' :

        # THESE MODES ARE STILL UNDER DEVELOPMENT.

        detector, aperture = _get_stis_imag_parameters(form_input)
        obsmode = make_stis_imaging_obsmode(detector, aperture)

        idict = dict(obsmode=obsmode,
                     detector=detector,
                     aperture=aperture,
                     crsplit=crsplit,
                     gain=gain,
                     binx=binx)

        idict["wodict"] = wodict = {}
        wodict["detector_type"] =  DETECTOR_TYPES[detector].lower()
        #idict['wodict'] = dict(detector_type=DETECTOR_TYPES[detector].lower())
        wodict["aperture"] = obsmode.split(',')[-1]
        idict["filter"] = wodict["aperture"].lower()
        wodict["obsmode"] = obsmode

    elif science_mode == 'targetacquisition':

        detector, aperture = _get_stis_imagacq_parameters(form_input)
        #the default extraction region has to be set to 32x32 for a PEAK ACQ aperture
        #this is currently being done by a hack in the extraction.py code which needs to
        #be changed at some point.

        aperture=form_input[aperture]
        ccdmode=form_input['ccdMode'].lower()
        point_target=(form_input['fsourceType'] == 'point')

        #The obsmode for imaging TA must specify the science aperture.
        obsmode = make_stis_targetacquisition_obsmode(detector,aperture,ccdmode,point_target)

        idict = dict(obsmode=obsmode,
                     detector=detector,
                     aperture=aperture,
                     ccdmode=ccdmode,
                     crsplit=crsplit,
                     gain=gain,
                     binx=binx)

        idict["wodict"] = wodict = {}
        wodict["detector_type"] = DETECTOR_TYPES[detector].lower()
        wodict["aperture"] = obsmode.split(',')[-1]
        idict["filter"] = wodict["aperture"].lower()
        wodict["obsmode"] = obsmode
        wodict["ccdmode"] = ccdmode

    else:
        raise ValueError('%s mode not implemented for stis'%science_mode)

    # this applies only to imaging and spectroscopic modes.
    idict["fuvglowregion"] = form_input.get("fuvglowregion", None)
    idict["ccddarklevel"] = form_input.get("ccddarklevel", None)

    return idict


def get_default(fdict, name):
    # Gets a default of 1 if name is not in the dictionary.
    try:
        return fdict[name]
    except KeyError:
        return 1


def make_stis_spectroscopic_obsmode(idict):
    disperser, detector, grating, central_wavelength, aperture = _get_stis_spec_parameters(idict)
            
    if central_wavelength != None:
        ci = CI_TABLE[grating][central_wavelength]
        obsmode = ",".join(['stis', detector, grating, ci + str(central_wavelength), aperture])

    else: #use default central wavelength if there is one
        try:
            grating_info = CI_TABLE[grating]
            central_wavelength = list(grating_info.keys())[0] #assume only 1 central wavelength in list
            designation=grating_info[central_wavelength]
            obsmode = ",".join(['stis', detector, grating, designation + str(central_wavelength), aperture])
        except KeyError:
            #otherise leave it off
            obsmode = ",".join(['stis', detector, grating, aperture])
    return obsmode


APERTURE_KEY = {
    "g750l" : "ccdaperture0",
    "g750m" : "ccdaperture0",
    "g430l" : "ccdaperture0",
    "g430m" : "ccdaperture0",
    "g705l" : "ccdaperture0",
    "g230lb" : "ccdaperture0",
    "g230mb" : "ccdaperture0",
    "g140l" : "fuvmamaaperture0",
    "g140m" : "fuvmamaaperture0",
    "g230l" : "nuvmamaaperture0",
    "g230m" : "nuvmamaaperture0",
    "e230m" : "e230maperture0",
    "e230h" : "e230haperture0",
    "e140m" : "e140maperture0",
    "e140h" : "e140haperture0",
    "prism" : "prismaperture0",
    }

def _get_stis_spec_parameters(idict):
    disperser = idict['disperser']
    detector, grating = disperser.split('_0',1)
    aperture = idict[APERTURE_KEY[grating]]
    try:
        central_wavelength = idict[grating.upper() + '_CentralWavelength']
    except KeyError:
        central_wavelength = None
    return disperser, detector, grating, central_wavelength, aperture

def make_stis_imaging_obsmode(detector, aperture=None):
    items = ['stis',detector,'mirror']
    if aperture is not None:
        items.append(aperture)
    obsmode = ",".join(items)
    return obsmode

def make_stis_targetacquisition_obsmode(detector,aperture,ccdmode,point_target):

    if "x" in aperture.lower() and "peak" in ccdmode.lower():
        #we need to put an s and remove the decimals
        dloc=aperture.split(".")
        aperture='s'+''.join(dloc)
    
    if "peak" in ccdmode.lower():
        if point_target:
            items = ['stis',detector,aperture]
        else:
            if "nd" in aperture.lower():
                items = ['stis',detector,aperture]
            else:
                items = ['stis',detector]        
    else:
        items = ['stis',detector,aperture]
        
    obsmode = ",".join(items)

    return obsmode

def _get_stis_imag_parameters(idict):
    detector = idict['detector']
    aperture = idict[detector + 'aperture0']
    return detector, aperture

def _get_stis_imagacq_parameters(idict):
    if idict["ccdMode"] == "ACQ":
        aperture = "ccdapertureACQ"
    else:
        aperture = "ccdapertureACQPEAK"
    detector = idict["detector"] = "ccd"
    return detector, aperture
