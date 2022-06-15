#TODO web cleanup:
# there is some duplication in between WFC3 UVIS and IR modules.


from __future__ import division

import pandeia.engine.helpers.bit.instruments
import pandeia.engine.helpers.bit.pyetc_util as util

#
# read a dictionary from one of the .dat files
def get_data( name ) :
    return pandeia.engine.helpers.bit.instruments.get_instrument_data( 'wfc3uvis', name )

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

    idict = dict( detector=form_input['detector'].lower(),
                  nreads = form_input['crsplit'],
                  )
    #check that crsplit is valid
    if idict['nreads'] < 1:
        raise ValueError("CRSPLIT must be greater than 0")

    #Web output dictionary that will be stuffed into db
    idict['wodict'] = wodict = util.stable_dict()
    #all wfc3uvis detectors are ccd-like
    wodict['detector_type'] = 'ccd'

    # Only show post flash in certain situations, see #922
    wodict["show_post_flash"] = False

    if science_mode == 'imaging' or science_mode == 'scimaging':
        #make an obsmode string that can be submitted to pysynphot
        idict['obsmode'], myfilter = make_wfc3uvis_imaging_obsmode(form_input)
        wodict['filter_info']=FILTERS[myfilter]
        wodict['filter_name']=myfilter
        idict['filter'] = myfilter.lower()
        idict['post_flash_electrons'] = form_input['post_flash_wfc3']
        wodict["show_post_flash"] = True

        return idict

    if science_mode == 'spectroscopic':
        #make an obsmode string that can be submitted to pysynphot
        idict['obsmode'] = make_wfc3uvis_spectroscopic_obsmode(form_input)
        idict['obswave']=form_input['obswave']
        idict['filter']='g280'
        idict['grating']='g280'
        return idict
    else:
        raise ValueError('%s mode not implemented for wfc3uvis'%science_mode)
    return idict


def make_wfc3uvis_imaging_obsmode(form_input):

    # figure out what kind of filter was selected
    ftype = form_input['wfc3_filter_type']
    myfilter = form_input[FILTER_MAP[ftype]]
    obsmode='wfc3,'+form_input['detector'].lower() + ','  + myfilter
    return obsmode, myfilter

def make_wfc3uvis_spectroscopic_obsmode(form_input):
    return  'wfc3,'+form_input['detector'].lower() + ',g280'
