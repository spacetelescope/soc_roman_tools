from __future__ import division

import pandeia.engine.helpers.bit.instruments
#
# read a dictionary from one of the .dat files
def get_data( name ) :
    return pandeia.engine.helpers.bit.instruments.get_instrument_data( 'cos', name )

CONFIG = get_data("config")

def instrument_from_form(science_mode, form_input):
    """Function to extract values of interest from the web form
    data and return a smaller dictionary which will be used to
    update the engine_input dictionary.
    """
    if science_mode == 'spectroscopic' or science_mode == 'spectroscopicacq':
        obsmode = make_cos_spectroscopic_obsmode(form_input)
        disperser, detector, grating, central_wavelength, aperture = _get_cos_spec_parameters(form_input)
        idict = dict(disperser=disperser,
                     detector=detector,
                     grating=grating,
                     central_wavelength=central_wavelength,
                     aperture=aperture,
                     obsmode=obsmode)
        idict["wodict"] = wodict = {}
        # all cos detectors are mamas
        wodict["detector_type"] = 'mama'

        wodict["grating_info"] = {}
        wodict["grating_info"]["resolution"] = CONFIG[detector]["grating_resolution"][grating]
        
        if science_mode == 'spectroscopic':
            idict['obswave'] = form_input['obswave']
        else:
            #enforce a specacq prohibition
            if (form_input['disperser'] == 'nuv_0g230l' and
                form_input['nuvMode'] == 'ACQ/PEAKXD' and
                form_input['Stripe'] == 'C'):
                raise ValueError("You have selected grating G230L with Stripe C. This is not a valid mode for Spectroscopic Target Acquisition/PeakXD.")
            
    elif science_mode == 'imaging' or science_mode == 'targetacquisition':
        obsmode = make_cos_imaging_obsmode(form_input)
        mirror, aperture = _get_cos_imag_parameters(form_input)
        detector = 'nuv'
        idict = dict(obsmode=obsmode,
                     detector=detector,
                     mirror=mirror,
                     aperture=aperture)
        idict["wodict"] = wodict = {}
        # all cos detectors are mamas
        wodict["detector_type"] = 'mama'
        wodict["mirror"] = {"mirrora" : "Mirror A", "mirrorb" : "Mirror B"}[mirror]
        idict["filter"] = mirror
    else:
        raise ValueError('%s mode not implemented for cos'%science_mode)
    wodict["aperture"] = {
        "psa" : "Primary Science Aperture",
        "boa" : "Bright Object Aperture"
    }[aperture]
    return idict

def make_cos_spectroscopic_obsmode(idict):
    disperser, detector, grating,  central_wavelength, aperture = _get_cos_spec_parameters(idict)
    obsmode = ",".join(['cos', detector, grating, 'c'+str(central_wavelength), aperture])
    return obsmode

def make_cos_spectroscopicacq_obsmode(idict):
    return make_cos_spectroscopic_obsmode(idict)

def _get_cos_spec_parameters(idict):
    disperser = idict['disperser']
    detector, grating = disperser.split('_0',1)
    central_wavelength = idict[grating.upper() + '_CentralWavelength']
    aperture = idict['cosaperture0'].lower()
    return disperser, detector, grating, central_wavelength, aperture

def make_cos_imaging_obsmode(idict):
    mirror, aperture = _get_cos_imag_parameters(idict)
    obsmode = ",".join(['cos',mirror,'nuv',aperture])
    return obsmode

def make_cos_targetacquisition_obsmode(idict):
    return make_cos_imaging_obsmode(idict)

def _get_cos_imag_parameters(idict):
    mirror = idict['Mirror']
    aperture = idict['cosaperture0'].lower()
    return mirror, aperture
