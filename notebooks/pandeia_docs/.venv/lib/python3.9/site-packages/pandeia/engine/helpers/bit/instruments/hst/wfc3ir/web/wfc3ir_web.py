#TODO web cleanup:
# there is some duplication in between WFC3 UVIS and IR modules.


from __future__ import division

import pandeia.engine.helpers.bit.pyetc_util as util
import pandeia.engine.helpers.bit.instruments

#
# read a dictionary from one of the .dat files
def get_data( name ) :
    return pandeia.engine.helpers.bit.instruments.get_instrument_data( 'wfc3ir', name )

FILTERS = get_data("config")["filters"]

FILTER_MAP = {
    'narrow': 'wfc3_filter_n',
    'medium': 'wfc3_filter_m',
    'wide':   'wfc3_filter_w',
    'quad':   'wfc3_filter_q',
}


def instrument_from_form(science_mode, form_input):
    """Function to extract values of interest from the web form
    data and return a smaller dictionary which will be used to
    update the engine_input dictionary.
    """

    idict = dict( nreads = form_input['crsplit'],
                  detector='ir'
                  )
    #check that crsplit is valid
    if idict['nreads'] < 1:
        raise ValueError("CRSPLIT must be greater than 0")

    #Web output dictionary that will be stuffed into db
    idict['wodict'] = wodict = util.stable_dict()
    #the wfc3 IR detectors are ccd-like
    wodict['detector_type'] = 'ccd'
    wodict['IR'] = True
    wodict['HeliumStandard'] = form_input['HeliumStandard']

    if science_mode == 'imaging' or science_mode == 'scimaging':
        #make an obsmode string that can be submitted to pysynphot
        idict['obsmode'],myfilter =make_wfc3ir_imaging_obsmode(form_input)
        wodict['filter_info']=FILTERS[myfilter]
        wodict['filter_name']=myfilter
        idict['filter'] = myfilter

        return idict

    if science_mode == 'spectroscopic' or science_mode == 'scspectroscopic':
        #convert the obswave in the form from microns to angstroms
        idict['obswave']=form_input['obswave'] * 10000.
        idict['grating']=form_input['disperser'].lower()[-4:]

       #make an obsmode string that can be submitted to pysynphot
        idict['obsmode'] = make_wfc3ir_spectroscopic_obsmode(form_input)

        return idict
    else:
        raise NotImplementedError('%s mode not implemented for wfc3ir'%science_mode)
    return idict


def make_wfc3ir_imaging_obsmode(form_input):
    #Based on the detector, we decide which filters to pick out.
    myfilter = form_input['irfilt0'].lower()
    obsmode='wfc3,ir,'+myfilter
    return  obsmode, myfilter

def make_wfc3ir_spectroscopic_obsmode(form_input):
    myfilter = form_input['disperser'].lower()[4:]

    #this mode ONLY has grisms, so assume filter is a grism to start with
    obsmode='wfc3,ir,'  + myfilter
    return obsmode
