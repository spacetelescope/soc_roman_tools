#!/usr/bin/env python

# basics
import glob, json, os, pprint, sys
import numpy as np

# pyetc tools  (hopefully not too many here)
import pysynphot
from pysynphot import spparser as SP
import synphot as syn
from astropy import units as u
from astropy import constants as const

# pandeia tools
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation
from pandeia.engine.io_utils import NumPyArangeEncoder
from pandeia.engine.instrument_factory import InstrumentFactory
from pandeia_test import utils # for rglob - maybe eventually move that to pandeia


# doc
"""
The general idea here is to take in an HST ETC (HETC) engine test file i.e. a
*pyetc* engine input file (*.peng), and convert it to *pandeia* engine inputs
outputting a '*.jeng' file.

We are starting with:
- HST inputs → JWST inputs
- JWST outputs → HST outputs (later)
- assumptions at start:
  - turn the sky off
  - simmode = SNR
"""

DATA_DIR = os.environ['pandeia_refdata']

# data drive this whole thing
HST_INSTS = ['acs','cos','stis','wfc3ir','wfc3uvis']
HST_INSTS_NEW_NAMES = ['wfc3']
HST_DETECTORS = ['ccd','ir','fuv','nuv','fuvmama','nuvmama','sbc','uvis1','uvis2','wfc', 'hrc']
HST_MODES = ['imaging','spectroscopic', 'spectroscopicacq', 'targetacquisition', 'scimaging', 'scspectroscopic', 'rampfilter']
HST_MIRRORS = ['mirrora','mirrorb']

MODE_NAME_DICT = {'spectroscopic': 'spectroscopy', 'scspectroscopic': 'scan_spectroscopy', 'scimaging': 'scan_imaging'}

DETECTOR_INSTRUMENT_DICT = {'wfc': 'acs',
                            'hrc': 'acs',
                            'sbc': 'acs',
                            'fuv': 'cos',
                            'nuv': 'cos',
                            'ccd': 'stis',
                            'fuvmama': 'stis',
                            'nuvmama': 'stis',
                            'ir': 'wfc3',
                            'uvis1': 'wfc3',
                            'uvis2': 'wfc3'}


# Dictionaries to convert the listed flux/wavelength units in the peng file to their corresponding astropy/pysynphot units.
FLUX_UNIT_DICT = {'photlam':syn.units.PHOTLAM,
                  'photnu':syn.units.PHOTNU,
                  'flam':syn.units.FLAM,
                  'fnu':syn.units.FNU,
                  'jy':u.Jy,
                  'stmag':u.STmag,
                  'abmag':u.ABmag}

WAVE_UNIT_DICT = {'angstroms':u.angstrom,
                  'nm':u.nm,
                  'm':u.m,
                  'micron':u.micron,
                  'microns':u.micron
                  }

SOURCEVAL_FUNC_DICT = {"spec":"sed",
                       "z":"redshift",
                       "icat":"sed",
                       "unit":"sed",
                       "pl":"sed",
                       "bb":"sed",
                       "em":"lines",
                       "rn":"normalization",
                       "ebmvx":"extinction"}

PWEB_FILTER_DICT = {"Johnson/V":"johnson,v", 
                    "Johnson/U":"johnson,u",
                    "Johnson/B":"johnson,b", 
                    "Johnson/R":"johnson,r", 
                    "Johnson/I":"johnson,i",
                    "Johnson/J":"johnson,j",
                    "Johnson/H":"johnson,h",
                    "Johnson/K":"johnson,k",
                    "Cousins/R":"cousins,r",
                    "Cousins/I":"cousins,i", 
                    "Bessell/J":"bessell,j", 
                    "Bessell/H":"bessell,h", 
                    "Bessell/K":"bessell,k",
                    "Sloan/U":"sdss,u",
                    "Sloan/G":"sdss,u",
                    "Sloan/R":"sdss,u",
                    "Sloan/I":"sdss,u",
                    "Sloan/Z":"sdss,u",
                    "Galex/FUV":"galex,fuv",
                    "Galex/NUV":"galex,nuv",
                    "NICMOS/F110W":"nicmos,2,f110w",
                    "NICMOS/F160W":"nicmos,2,f160w",
                    "ACS/F435W":"acs,wfc1,f435w",
                    "ACS/F475W":"acs,wfc1,f475w",
                    "ACS/F555W":"acs,wfc1,f555w",
                    "ACS/F606W":"acs,wfc1,f606w",
                    "ACS/F625W":"acs,wfc1,f625w",
                    "ACS/F775W":"acs,wfc1,f775w",
                    "ACS/F814W":"acs,wfc1,f814w",
                    "ACS/F850LP":"acs,wfc1,f850lp",
                    "WFC3/UVIS/F218W":"wfc3,uvis1,f218w",
                    "WFC3/UVIS/F200LP":"wfc3,uvis1,f200lp",
                    "WFC3/UVIS/F225W":"wfc3,uvis1,f225w",
                    "WFC3/UVIS/F275W":"wfc3,uvis1,f275w",
                    "WFC3/UVIS/F300X":"wfc3,uvis1,f300x",
                    "WFC3/UVIS/F336W":"wfc3,uvis1,f336w",
                    "WFC3/UVIS/F350LP":"wfc3,uvis1,f350lp",
                    "WFC3/UVIS/F390W":"wfc3,uvis1,f390w",
                    "WFC3/UVIS/F438W":"wfc3,uvis1,f438w",
                    "WFC3/UVIS/F475X":"wfc3,uvis1,f475x",
                    "WFC3/UVIS/F600LP":"wfc3,uvis1,f600lp",
                    "WFC3/IR/F098M": "wfc3,ir,f098m",
                    "WFC3/IR/F105W": "wfc3,ir,f105w",
                    "WFC3/IR/F110W": "wfc3,ir,f110w",
                    "WFC3/IR/F125W": "wfc3,ir,f125w",
                    "WFC3/IR/F140W": "wfc3,ir,f140w",
                    "WFC3/IR/F160W": "wfc3,ir,f160w"}


PWEB_NORMUNIT_DICT = {"angstroms": "angstrom",
                      "nanometers": "nm",
                      "microns": "micron"
                        }

# unlike peng2jeng, we are not tying ourselves to only translate and produce currently-valid ETC calculations.
# though we should have a flag that will only print valid calculations as determined by successful translation 
# and/or appearance of the value in the relevant config.json

# In this new construction, the second item in the list is either the location in the .jeng dictionary OR 
# the name of a function that will parse that input.
# The same construction as before holds: 
# 1. We load the .peng file, then parse each line to build up a list of what to parse.
# 2. Then we create a default .jeng calculation dictionary
# 3. All the simple translations are just added to the dictionary
# 4. We go through all the '_parse' functions to populate things that need translation or to read from multiple fields
# 5. Then we go through the '_valid' functions to check the original value against the new one.

def get_instdotmode_from_output(str_blob):
    ''' grab the instrument+mode value out of the output - bit of a tiny hack for reporting '''
    assert 'Wrote ' in str_blob
    assert ' test file: ' in str_blob
    for line in str_blob.split('\n'):
        if line.startswith('Wrote ') and ' test file: ' in line:
            # grab idm = inst-dot-mode (e.g. "stis.echelle")
            idm = line.split(' ')[1]
            inst = idm.split('.')[0]
            #assert inst in HST_INSTS + HST_INSTS_NEW_NAMES, 'unexpected inst name: %s' % inst
            return idm
    assert False, 'Expected string not found in the str_blob: %s' % str_blob


def safe_eval(astr):
    '''
    this is a stand-in for pyetc.main.util:safe_eval() because we do NOT want to pull all of
    that in here simply to load a dictionary file...
    '''
    lines = astr.split('\n')
    nocomments = [l for l in lines if not l.strip().startswith('#')]
    buf = ' '.join(nocomments)
    buf = buf.replace('\t',' ') # the result is valid Python, but not yet valid JSON

    # NOTE - since this is not public facing, we are going to use eval here.  This
    # is a temporary tool and will be dropped after the conversion to pandeia.
    return eval(buf)


def _safe_update(insdict, key, value):
    if key in insdict:
        assert insdict[key] == value, f"Inconsistent values for {key}: {insdict[key]}, {value} "
    insdict[key] = value
    return insdict

def pyetc_main_util_read_dict(fname):
    '''
    (stolen directly from pyetc.main.util so as to not require a depenency on pyetc if possible)
    read a python dictionary from a file that was written with write_dict.
    '''
    # no six for this
    if sys.version_info[0] >= 3:
        f=open(fname,'r', encoding="utf-8")
    else:
        f=open(fname,'r')
    datastr = f.read()
    f.close()
    # convert DOS file to Unix - otherwise the eval will fail
    datastr = datastr.replace('\r','')
    try :
        datadict = safe_eval(datastr)
    except Exception as e:
        raise ValueError('Exception: %s, from data file: %s' % (e, fname))
    return datadict


def peng2dict(fname):
    " convenience function to load an input file "
    assert os.path.exists(fname), 'Input file unfound: %s' % fname
    # can't just use json.load for pyetc peng files
    return pyetc_main_util_read_dict(fname)

# The correctors, for complex translations

def disperser_valid(adict, calculation):
    assert calculation['configuration']['instrument']['disperser'] in adict['disperser'] or calculation['configuration']['instrument']['aperture'] in adict['disperser'], f"Inconsistent disperser config {adict['disperser']}"
    return calculation

def detector_parse(adict, calculation, pweb_key):
    # We need to settle once and for all if we're making detector independent of aperture. This and obsmode_parse disagree on that point.
    # Basically, we are... except for WFC3 UVIS1/UVIS2, which are so similar they will share modes.

    # See if we need to translate the mode name, otherwise use what's there already.
    mode_name = MODE_NAME_DICT.get(adict['science_mode'], adict['science_mode'])

    if 'uvis' in adict['extraction']['detector']:
        calculation['configuration']['instrument']['mode'] = f"uvis_{mode_name}"
    else:
        calculation['configuration']['instrument']['mode'] = f"{adict['extraction']['detector']}_{mode_name}"


    calculation['configuration']['instrument']['aperture'] = adict['extraction']['detector']
    if 'aperture' in adict:
        if adict['aperture'] != adict['extraction']['detector'] and adict['instrument'] != 'stis':
            # STIS pengs have the filter in the aperture keyword.  Just leave aperture as extraction/detector.
            calculation['configuration']['instrument']['aperture'] = adict['aperture']
    assert adict['detector'] == adict['extraction']['detector']
    return calculation

def obswave_parse(adict, calculation, pweb_key):
    # Convert the value from Angstroms to microns and store in ['strategy']['reference_wavelength']
    val = float(adict['obswave']) * u.angstrom
    calculation['strategy']['reference_wavelength'] = val.to_value(u.micron)

    return calculation

def obsmode_parse(adict, calculation, pweb_key):

    # Instrument-dependent
    # science mode dependent
    # 
    pweb_dict = locate_pweb(pweb_key)
    obsmode_list = adict['obsmode'].split(',')
    inst = calculation['configuration']['instrument']['instrument']
    if 'wfc3' in inst:
        inst = 'wfc3'

    assert inst in obsmode_list[0], "Invalid instrument"
    result = False
    # some STIS obsmode strings have neither filter nor disperser
    # for item in obsmode_list:
    #     if calculation['configuration']['instrument']['filter'] is not None and calculation['configuration']['instrument']['filter'] in item.lower():
    #         result = True
    #     if calculation['configuration']['instrument']['disperser'] is not None and calculation['configuration']['instrument']['disperser'] in item.lower():
    #         result = True
    # assert result

    if "MJD" in obsmode_list[-1]:
        """
        MJD is not required, but is present in virtually all pyetc tests (added 
        by pyetc/etc_web/etc/perform_calculation.py _handle_mjd(), called by
        process_form_input().
        All further notes will consider this to have been removed.
        """
        calculation['configuration']['instrument']['mjd'] = int(obsmode_list.pop(-1).split("#")[1])
    if 'cos' in calculation['configuration']['instrument']['instrument']:
        if 'spectro' in adict['science_mode']:
            """
            instrument,detector,grating,cenwave,aperture
            Identical processing to stis spectroscopy except that the cenwave is always 
            present, and prepended with ‘c’. There is some hardcoded checking in the 
            instrument_from_form function for invalid combinations of grating and stripe
            """
            calculation['configuration']['instrument']['detector'] = obsmode_list[1].lower()
            calculation['configuration']['instrument']['disperser'] = obsmode_list[2].lower()
            calculation['configuration']['instrument']['cenwave'] = obsmode_list[3].lower()
            calculation['configuration']['instrument']['aperture'] = obsmode_list[4].lower()
        elif ('targetacquisition' in adict['science_mode']) or ('imaging' in adict['science_mode']):
            """
            instrument,mirror,detector,aperture
            set by make_cos_imaging_obsmode, all of the heavy lifting is done by 
            _get_cos_imag_parameters, which pulls mirror out of 'Mirror' in the pweb file, 
            and aperture out of ‘cosaperture0’ in the pweb file. detector is present and 
            always ‘nuv’. There are no filters.
            """
            calculation['configuration']['instrument']['detector'] = obsmode_list[2].lower()
            calculation['configuration']['instrument']['aperture'] = obsmode_list[1].lower()
            calculation['configuration']['instrument']['slit'] = obsmode_list[3].lower()
            calculation['configuration']['instrument']['filter'] = None
            calculation['configuration']['detector']['nsplit'] = 1
        else:
            raise NameError(f"Could not identify science_mode {adict['science_mode']}")
    elif 'stis' in calculation['configuration']['instrument']['instrument']:
        if 'spectro' in adict['science_mode']:
            """
            instrument,detector,disperser,[cenwave],aperture
            All work handled by make_stis_spectroscopic_obsmode, which calls 
            _get_stis_spec_parameters. The thing that calls make_stis_spectroscopic_parameters 
            then reruns _get_stis_spec_parameters to get values for itself. But! We then pull 
            the aperture values out of the OBSMODE to put in wo_dict; it’s not modified by 
            make_stis_spectroscopic_obsmode so I don’t know why.
            The only optional component is cenwave, which make_stis_spectroscopic_obsmode 
            modifies by prepending a c or i depending on a hardcoded table; cenwave is either 
            [3] or not present. Aperture is always last.
            """
            calculation['configuration']['instrument']['detector'] = obsmode_list[1].lower()
            calculation['configuration']['instrument']['disperser'] = obsmode_list[2].lower()
            calculation['configuration']['instrument']['aperture'] = obsmode_list[-1].lower()
            if len(obsmode_list) > 4:
                calculation['configuration']['instrument']['cenwave'] = obsmode_list[3].lower()
        elif 'targetacquisition' in adict['science_mode']:
            """
            stis targetacquisition: Basically same as imaging, except there is no "mirror" 
            item. And if the ccdmode is peak and the aperture/filter has an x in it, we 
            remove any '.' in the aperture string and prepend an s. Also, there is no 
            “mirror” item. When acq is peak, the aperture IS an aperture (or rather an 
            aperture/filter combo, given that they seem to have ND filters); when acq 
            isn’t peak, it’s a filter.
            """
            calculation['configuration']['instrument']['detector'] = obsmode_list[1].lower()
            if len(obsmode_list) > 2:
                if pweb_dict:
                    if 'peak' in pweb_dict['ccdMode'].lower():
                        calculation['configuration']['instrument']['aperture'] = obsmode_list[2].lower()
                    else:
                        calculation['configuration']['instrument']['filter'] = obsmode_list[2].lower()
                else:
                    raise FileNotFoundError("Could not find matching pweb file")
        elif 'imaging' in adict['science_mode']:
            """
            instrument,detector,”mirror”,[aperture]
            The work is done by _get_stis_imag_parameters, which determines detector and 
            aperture (which is apparently actually the filter in idict; it’s the aperture 
            in wodict). The other parameter (possibly the actual aperture? or grating? is 
            “mirror”, hardcoded. _stis_imag_obsmode implies that aperture (aka filter) is 
            optional, but given that even the 50CCD/clear position is specified in PyETC 
            obsmode strings (and you can't use an ACTUAL slit/aperture), I suspect it’s 
            NOT optional (for PyETC use, anyway)
            """
            calculation['configuration']['instrument']['detector'] = obsmode_list[1].lower()
            # STIS pengs have the filter in the aperture and obsmode settings.
            #  Aperture is set correctly in the detector parsing function.
            if len(obsmode_list) > 3:
                calculation['configuration']['instrument']['filter'] = obsmode_list[3].lower()
        else:
            raise NameError(f"Could not identify science_mode {adict['science_mode']}")
    elif 'wfc3' in calculation['configuration']['instrument']['instrument']:
        if 'spectro' in adict['science_mode']:
            """
            wfc3ir
            instrument,ir,disperser
            disperser (called filter in the code) comes from the pweb file as 
            ‘disperser’ (made lowercase). The other two items are hardcoded. The 
            only complexity is that we strip the first four characters from 
            ‘disperser’ in the form to get the disperser. The other oddity is that 
            the “grating” in the peng file is the LAST four characters of the 
            pweb[‘disperser’] entry. A quick check of a pweb file shows that the 
            disperser contains an 8-char string (ir_0G102 or similar) so both bits 
            of code just happen to be doing the same thing.
            wfc3uvis
            instrument,detector,g280
            Done by make_wfc3uvis_spectroscopic_obsmode; it only pulls ‘detector’ 
            from the form_input and hardcodes the g280 grating. (it IS the only 
            one, and not selectable)
            """

            calculation['configuration']['instrument']['detector'] = obsmode_list[1].lower()
            calculation['configuration']['instrument']['disperser'] = obsmode_list[2].lower()
        elif "imaging" in adict['science_mode']:
            """
            wfc3ir
            instrument,ir,filter
            filter comes directly from the pweb file as ‘irfilt0’ (made lowercase). 
            The other two items are hardcoded. It’s that simple.
            wfc3uvis
            instrument,detector,filter
            work is done by make_wfc3uvis_imaging_obsmode, which checks the 
            wfc3_filter_type in the form to figure out which form input to read for 
            the filter. Detector is ‘detector’ from the form input.
            """
            calculation['configuration']['instrument']['detector'] = obsmode_list[1].lower()
            calculation['configuration']['instrument']['filter'] = obsmode_list[2].lower()
        else:
            raise NameError(f"Could not identify science_mode {adict['science_mode']}")
    elif 'acs' in calculation['configuration']['instrument']['instrument']:
        if 'spectro' in adict['science_mode']:
            """
            instrument,detector,[disperser],[filter]
            Complex selection methods, because this is sometimes called with an input 
            (peng) dictionary instead of a pweb dictionary; sometimes the detector and 
            grating are separate items (peng), sometimes they’re combined into a single 
            ‘disperser’ string (pweb) with detector_disperser, and the form has to be 
            parsed to figure out which detector’s filter to pull out (with the 
            additional complication that hrc’s listings are backwards - if the disperser 
            starts with 0, look at hrcfilt1; if the disperser starts with 1, look at 
            hrcfilt0). There can either be both or neither of them if any are “clear”. 
            (though the web interface doesn’t allow the disperser to be clear, so it’s 
            effectively required.) Left unsaid: some of the filters are actually 
            polarizers. (of course, given that neither ETC intends to actually implement 
            polarization, they are effectively filters)
            """
            calculation['configuration']['instrument']['detector'] = obsmode_list[1].lower()
            if obsmode_list[1] == 'wfc1':
                calculation['configuration']['instrument']['detector'] = 'wfc'

            # According to ACS Handbook, dispersers start with 'G8' or 'PR'.
            # Polarizers start with 'POL'
            # Everything else will be assumed to be a filter.

            # Loop through the remaining obsmode_list items, and determine if they are filter/polarizer/disperser...
            for i in range(len(obsmode_list) - 2): #[0] and [1] are instrument and detector, so skip the first 2 items.
                index = i + 2
                if obsmode_list[index][:2] in ['g8', 'pr']:
                    calculation['configuration']['instrument']['disperser'] = obsmode_list[index].lower()
                elif str.startswith(obsmode_list[index], 'pol'):
                    calculation['configuration']['instrument']['polarizer'] = obsmode_list[index].lower()
                else:
                    calculation['configuration']['instrument']['filter'] = obsmode_list[index].lower()
        elif 'ramp' in adict['science_mode']:
            """
            instrument,detector,filter#wave
            As with acs imaging, the ETC assumes any use of wfc is using wfc1 as the detector. 
            Filter is read out of the pweb file as ‘(detector)filt1’, and the obswave is read 
            from ‘obswave’ and concatenated. The problem here is that pysynphot has a bug (at 
            least as of 7 years ago) in that it couldn’t support the full wavelength range of 
            the ramp filters, and the solution was to actually run pysynphot on the obsmode 
            string to see if it crashed, and then tell the user the wavelength is out of 
            range. Presumably this script will never run into one of those because those .peng
            files could never be made.
            """
            calculation['configuration']['instrument']['detector'] = obsmode_list[1].lower()
            if obsmode_list[1] == 'wfc1':
                calculation['configuration']['instrument']['detector'] = 'wfc'
            calculation['configuration']['instrument']['filter'] = obsmode_list[2].split('#')[0].lower()
            calculation['strategy']['reference_wavelength'] = (float(obsmode_list[2].split('#')[1]) * u.angstrom).to_value(u.micron)
        elif 'imaging' in adict['science_mode']:
            """
            instrument,detector,[filter1],[filter2]
            Very complex processing, mostly done in make_acs_imaging_obsmode. wfc is split 
            across two detectors, but ETC assumes any wfc is using the wfc1 detector. hrc 
            and wfc can have two filters (taken from ‘(detector)filt0’ and ‘(detector)filt1’ 
            in the pweb dictionary. Of the input filters, the next element is any non-clear 
            filter(s), though later code throws an error if you have more than one non-clear 
            filter (which already shouldn’t be in the list of filters you created earlier) 
            There does not seem to be a prescribed order of filter and polarizer or filter 
            and clear (though it’s likely that, because all the polarizers are on wheel 2, 
            that they’re all the second item); a separate block of code finds and returns 
            them to the acs instrument_from_form() function. That block of code does make 
            it clear that all polarizers start with “pol” and all clear filters start with 
            “clear”. But again, the code should be weeding out clear filters, so you could 
            have a case where both are clear and the obsmode string ends at “detector”. Using 
            clear in both wheels is allowed, but generates a not-supported warning.
            """
            calculation['configuration']['instrument']['detector'] = obsmode_list[1].lower()
            if calculation['configuration']['instrument']['detector'] == 'wfc1':
                calculation['configuration']['instrument']['detector'] = 'wfc'
            if len(obsmode_list) == 2:
                calculation['configuration']['instrument']['filter'] = adict['filter'].lower()
            if len(obsmode_list) == 3:
                if "pol" in obsmode_list[2]:
                    calculation['configuration']['instrument']['polarizer'] = obsmode_list[2].lower()
                else:
                    calculation['configuration']['instrument']['filter'] = obsmode_list[2].lower()
            if len(obsmode_list) == 4:
                if "pol" in obsmode_list[3]:
                    calculation['configuration']['instrument']['filter'] = obsmode_list[2].lower()
                    calculation['configuration']['instrument']['polarizer'] = obsmode_list[3].lower()
                else:
                    calculation['configuration']['instrument']['filter'] = obsmode_list[3].lower()
                    calculation['configuration']['instrument']['polarizer'] = obsmode_list[2].lower()

    return calculation


def extraction_parse(adict, calculation, pweb_key):
    calculation['strategy']['gyromode'] = adict['extraction']['gyromode']
    calculation['configuration']['instrument']['detector'] = adict['extraction']['detector']
    calculation['strategy']['background_subtraction'] = False
    del calculation['strategy']['sky_annulus']

    # Pandeia expects the aperture_size parameter to be in arcsec.
    #  We need to read in the plate scale to convert pixels to arcsec.
    detector = adict['extraction']['detector']
    instrument = DETECTOR_INSTRUMENT_DICT[detector]
    with open(DATA_DIR + '/hst/%s/config.json' % instrument, 'r') as config_file:
        config_data = config_file.read()

    config = json.loads(config_data)

    # To get the plate scale, we need to know the aperture.  To get the aperture, set up an
    #  Instrument object.
    inst_config = {}
    inst_config['instrument'] = {}
    inst_config['instrument']['instrument'] = instrument
    # WFC3 UVIS1 and UVIS2 are so similar, they will share a mode.
    # See if we need to translate the mode name, otherwise use what's there already.
    mode_name = MODE_NAME_DICT.get(adict['science_mode'], adict['science_mode'])

    if 'uvis' in adict['extraction']['detector']:
        inst_config['instrument']['mode'] = f"uvis_{mode_name}"
    elif instrument == 'stis' and 'spec' in adict['science_mode']:
        if adict['grating'].startswith('e'):
            # Echelle mode
            inst_config['instrument']['mode'] = f"{detector}_echelle"
        elif adict['grating'].startswith('g'):
            # 1st order spec mode.  Now determine if slit or slitless.
            slit = "slit"
            if adict['aperture'].startswith('f') or adict['aperture'] == "50ccd" or adict['aperture'] == "25mama":
                slit = "slitless"
            inst_config['instrument']['mode'] = f"{detector}_{mode_name}_{slit}"
    else:
        inst_config['instrument']['mode'] = f"{detector}_{mode_name}"

    if 'filter' in adict.keys():
        inst_config['instrument']['filter'] = adict['filter'].lower()

    # Creating a STIS FUVMAMA or CCD instrument requires an extra param.
    #  Give it a value here (doesn't matter which for conversion purposes).
    if 'detector' in adict.keys() and 'stis' in instrument:
        inst_config['detector'] = {}
        if 'fuv' in adict['detector']:
            inst_config['detector']['fuv_glow_region'] = 'low'
        elif 'ccd' in adict['detector']:
            inst_config['detector']['dark_level'] = 'low'

    inst = InstrumentFactory(config=inst_config)
    aperture = inst.get_aperture()

    plate_scale_x = config['aperture_config'][aperture]['plate_scale_x'] # arcsec / pixel
    plate_scale_y = config['aperture_config'][aperture]['plate_scale_y'] # arcsec / pixel

    # Use the shape parameter and the science_mode parameter to determine the resulting strategy in the JENG file.
    if adict['extraction']['shape'] == 'rectangle' and 'spectro' not in adict['science_mode']:
        calculation['strategy']['method'] = 'imagingrectphot'
        calculation['strategy']['aperture_size'] = [adict['extraction']['size'][0] * plate_scale_x, adict['extraction']['size'][1] * plate_scale_y]
        calculation['strategy']['display_string'] = "Imaging Rectangular Photometry"
    elif adict['extraction']['shape'] == 'rectangle':
        calculation['strategy']['method'] = 'specapphot' # note that SpecApPhot only accepts the half-height of the box, and relies on other means to get the length of the extraction box
        calculation['strategy']['aperture_size'] = (adict['extraction']['size'][0] * plate_scale_x, adict['extraction']['size'][1] * plate_scale_y)
        calculation['strategy']['display_string'] = "Aperture Spectral Extraction"
    if adict['extraction']['shape'] == 'square':
        calculation['strategy']['method'] = 'imagingrectphot'
        calculation['strategy']['aperture_size'] = (adict['extraction']['size'] * plate_scale_x,adict['extraction']['size'] * plate_scale_y)
        calculation['strategy']['display_string'] = "Imaging Rectangular Photometry"
    elif adict['extraction']['shape'] == 'circle':
        calculation['strategy']['method'] = 'imagingapphot'
        calculation['strategy']['aperture_size'] = float(adict['extraction']['size'])
    elif adict['extraction']['shape'] == 'percent':
        calculation['strategy']['method'] = 'imagingapphot'
        calculation['strategy']['aperture_size'] = float(adict['extraction']['size']) * plate_scale_x/10.
        # TODO: by default, the units are 'ee' (encircled energy) which is currently not implemented in pandeia.  Once those tickets are completed, we can get rid of the conversions here and just use 'ee'.
        calculation['strategy']['units'] = 'arcsec'
    if 'stripe' in adict['extraction']:
        calculation['strategy']['order'] = adict['extraction']['stripe']
    if adict['extraction']['stype'] == 'extended':
        # According to Ivo, sdiameter and fdiameter are the same, and pyetc uses compatibility code like this
        try:
            diameter = adict['extraction']['sdiameter']
        except KeyError:
            try:
                diameter = adict['extraction']['fdiameter']
            except KeyError:
                diameter = adict['fdiameter']
        calculation['scene'][0]['shape'] = {"geometry": "flat","major": diameter, "minor": diameter,
            "norm_method": "integ_infinity", "surf_area_units": None }

    return calculation

def simmode_valid(adict, calculation):
    # as per https://innerspace.stsci.edu/display/JEP/2021-11-04+Reverse+Calc+API+Meeting+notes

    calculation['configuration']['detector']['time'] = adict['time']
    calculation['configuration']['detector']['snr'] = adict['SNR']
    calculation['configuration']['detector']['calculate_snr'] = adict['simmode'] == 'SNR'

    return calculation

def locate_pweb(pweb_key):
    """
    Returns the pweb dict given a "pweb key" which is either a peng fname (for which it is assumed to
    have a matching pweb file that we can find), OR it is the pweb dict itself.  In the former
    case we find the file, read it and return the dict.  In the latter case we just return the dict.

    This is generalized to be multipurpose.

    File path searching assumes the peng is in pyetc/test/engine.

    Returns None (without error) if no such file is found.  Callers chack this.
    """
    if isinstance(pweb_key, dict):
        return pweb_key # this is the pweb dict

    pweb_fname = pweb_key.replace('test/engine', 'test/web').replace('.peng', '.pweb')
    pweb_dict = None
    try:
        # load the pweb file
        with open(pweb_fname, 'r') as infile:
            pweb_dict = eval(infile.read())
    except FileNotFoundError:
        pass # pweb_dict = None

    return pweb_dict

def sourceval_parse(adict, calculation, pweb_key):
    sky_expr = adict['source_expr']
    pweb_src_keys = ['spec(', 'z(', 'icat', 'unit(', 'pl(', 'bb(', 'em(', 'rn(']

    try:
        # If this source_expr contains anything that might be easier to parse with a pweb, try and locate the pweb.
        if any(key in sky_expr for key in pweb_src_keys):
            pweb_dict = locate_pweb(pweb_key)

            if 'spec(' in sky_expr:
                # A spectrum from CDBS
                # key seems to not be stored in the pweb.  Still need to parse this from the peng I guess.
                #  spec always seems to be the innermost statement in a nest, so this shouldn't be difficult.
                sedstart = sky_expr.find('spec(')
                sedend = sky_expr.find(')', sedstart)
                spec = sky_expr[sedstart:sedend].replace('spec(', '')
                parser = spec.split('$')
                if parser[0] == "crgridbz77":
                    spec = 'grid/bz77/' + parser[1]
                elif parser[0] == "crgrid":
                    spec = 'grid/' + parser[1]
                elif parser[0] == "crcalspec":
                    spec = 'calspec/' + parser[1]
                elif parser[0] == "": # starts with $PYSYN_CDBS
                    spec = "/".join(spec.split("/")[1:])
                calculation['scene'][0]['spectrum']['sed'] = {'sed_type': 'custom', 'key': spec}
            if 'z(' in sky_expr:
                calculation['scene'][0]['spectrum']['redshift'] = pweb_dict['fRedshift']
            if 'icat(' in sky_expr:
                # Comparing *189398.peng to *189398.pweb, it's not clear that, teff, metallicity, and log_g values are stored in the pweb.
                #  Parsing the old way for now.
                # A catalog grid to interpolate
                sedstart = sky_expr.find('icat(')
                sedend = sky_expr.find(')', sedstart)
                spec = sky_expr[sedstart:sedend].replace('icat(', '')
                spec_list = spec.split(',')
                calculation['scene'][0]['spectrum']['sed'] = {'sed_type': spec_list[0], 'teff': float(spec_list[1]),
                                                            'metallicity': float(spec_list[2]),
                                                            'log_g': float(spec_list[3])}
            if 'unit(' in sky_expr:
                calculation['scene'][0]['spectrum']['sed'] = {'sed_type': 'flat', 'unit': "flam" if pweb_dict['fIsLambda'] == "'true'" else "fnu"}
            if 'pl(' in sky_expr:
                calculation['scene'][0]['spectrum']['sed'] = {'sed_type': 'powerlaw', 'index': pweb_dict['fIndex'], 'unit': pweb_dict['rn_flux_lambda_units']}
            if 'bb(' in sky_expr:
                calculation['scene'][0]['spectrum']['sed'] = {'sed_type': 'blackbody', 'teff': pweb_dict['fbbtemp']}
            if 'em(' in sky_expr:
                # Multiple lines can be added.  Split by '+' and add a line for each 'em('
                #  Example expr: rn(unit(1.,flam),band(galex,fuv),25.000000,abmag)+em(2000,10,1.5e-13,flam)+em(2020,10,1.5e-13,flam)+em(2040,10,1.5e-13,flam)
                # Info not stored in pweb, parse it out the old way here.
                pieces = sky_expr.split('+')
                for piece in pieces:
                    if 'em(' in piece:
                        param_list = piece.replace('em(', '').replace(')', '').split(',')
                        line = {'center': (float(param_list[0])*u.angstrom).to_value(u.micron), 'emission_or_absorption': 'emission', 'profile': 'gaussian', \
                        'strength': syn.units.convert_flux(float(param_list[0]) * u.angstrom, float(param_list[2])*syn.units.FLAM, syn.units.FNU), 'width': (float(param_list[1])/2)/float(param_list[0]) * 3e8}
                        calculation['scene'][0]['spectrum']['lines'].append(line)
            if 'rn(' in sky_expr:
                calculation['scene'][0]['spectrum']['normalization']['norm_flux'] = float(pweb_dict['rn_flux_lambda'])
                calculation['scene'][0]['spectrum']['normalization']['norm_fluxunit'] = pweb_dict['rn_flux_lambda_units']

                if 'band' in sky_expr:
                    normtype = pweb_dict['fftype_filters'].split('.')[1]
                    calculation['scene'][0]['spectrum']['normalization']['bandpass'] = PWEB_FILTER_DICT[pweb_dict[f'filter.{normtype}']]
                    if 'acs' in normtype or "wfc3" in normtype or "nicmos" in normtype:
                        calculation['scene'][0]['spectrum']['normalization']['type'] = "hst"
                    else:
                        calculation['scene'][0]['spectrum']['normalization']['type'] = "photsys"
                    del calculation['scene'][0]['spectrum']['normalization']['norm_wave']
                    del calculation['scene'][0]['spectrum']['normalization']['norm_waveunit']
                    normstart = sky_expr.find('band(')
                if 'box' in sky_expr:
                    calculation['scene'][0]['spectrum']['normalization']['norm_wave'] = pweb_dict['rn_lambda'] #float(bandlist[0])
                    calculation['scene'][0]['spectrum']['normalization']['norm_waveunit'] = PWEB_NORMUNIT_DICT[pweb_dict['rn_lambda_units']] #'angstroms' or 'nanometers'
                    calculation['scene'][0]['spectrum']['normalization']['type'] = "at_lambda"
                    normstart = sky_expr.find('box(')
                if 'ebvmx' in sky_expr:
                    calculation['scene'][0]['spectrum']['normalization']['bandpass'] = PWEB_FILTER_DICT[pweb_dict['filter.ubvri']] #'v'
                    calculation['scene'][0]['spectrum']['extinction']['unit'] = pweb_dict['rn_flux_bandpass_units'] #'mag'
                    calculation['scene'][0]['spectrum']['extinction']['law'] = pweb_dict['febmvtype'] #"mwrv31"
                    calculation['scene'][0]['spectrum']['extinction']['value'] = pweb_dict['febv'] #val * 3.1

                    redstart = sky_expr.find('ebmvx(')
                    redend = sky_expr.find(')', redstart)

                    calculation['calculation']['settings'] = {}
                    if redend > normstart:
                        # I manually tested this, and it seems to work.  Changing the order of the reddening and band resulted in
                        #   different extinction_first values.
                        # if the reddening ended after the normalization band definition was started, we applied reddening later
                        # the normalization will have been removed when the reddening is defined, but this won't fail because the
                        # normalization entry isn't fully removed, and thus reddening couldn't start before it did.
                        calculation['calculation']['settings']['extinction_first'] = False
                    else:
                        calculation['calculation']['settings']['extinction_first'] = True
    except (FileNotFoundError, TypeError):
        calculation = old_sourceval_parse(adict, calculation) # the old parsing is inferior but does not require a PWEB file

    return calculation

def old_sourceval_parse(adict, calculation):
    sky_expr = adict['source_expr']

    sky_expr = sky_expr.split('+')
    for item in sky_expr:
        # basically, we're going to remove elements from this until we end up at just a renormalize statement
        if '*' in item:
            scaleindex = item.find('*')
            scalefactorstart = item.rfind('(', 0, scaleindex)
            scalefactor = item[scalefactorstart + 1:scaleindex].strip()

            # Figure out what function this scaling is being applied to.
            funcend = item.find('(', scaleindex)
            funcname = item[scaleindex+1:funcend].strip()
            calculation['scene'][0]['spectrum'][SOURCEVAL_FUNC_DICT[funcname]]['scale'] = scalefactor
            item = item.replace(scalefactor, '').replace('*', '')
        if 'spec(' in item:
            # A spectrum from CDBS
            sedstart = item.find('spec(')
            sedend = item.find(')', sedstart)
            spec_orig = item[sedstart:sedend+1]
            spec = item[sedstart:sedend].replace('spec(','')
            parser = spec.split('$')
            if parser[0] == "crgridbz77":
                spec = 'grid/bz77/' + parser[1]
            elif parser[0] == "crgrid":
                spec = 'grid/' + parser[1]
            elif parser[0] == "crcalspec":
                spec = 'calspec/' + parser[1]
            elif parser[0] == "": # starts with $PYSYN_CDBS
                spec = "/".join(spec.split("/")[1:])
            calculation['scene'][0]['spectrum']['sed'] = {'sed_type': 'cdbs', 'key': spec}
            item = item.replace(spec_orig,'')
        if 'z(' in item:
            zstart = item.find('z(')
            zend = item.find(')', zstart)
            z_orig = item[zstart:zend+1]
            z = item[zstart:zend].replace('z(', '').replace(',', '')
            calculation['scene'][0]['spectrum']['redshift'] = float(z)
            item = item.replace(z_orig, '')
        if 'icat(' in item:
            # A catalog grid to interpolate
            sedstart = item.find('icat(')
            sedend = item.find(')', sedstart)
            spec_orig = item[sedstart:sedend+1]
            spec = item[sedstart:sedend].replace('icat(','')
            speclist = spec.split(',')
            calculation['scene'][0]['spectrum']['sed'] = {'sed_type': speclist[0], 'teff': float(speclist[1]), 'metallicity': float(speclist[2]), 'log_g': float(speclist[3])}
            item = item.replace(spec_orig,'')
        if 'unit(' in item:
            # Flat spectrum
            sedstart = item.find('unit(')
            sedend = item.find(')', sedstart)
            spec_orig = item[sedstart:sedend+1]
            spec = item[sedstart:sedend].replace('unit(','')
            speclist = spec.split(',')
            if speclist[0] == '1.':
                calculation['scene'][0]['spectrum']['sed'] = {'sed_type': 'flat', 'unit': speclist[1]}
            item = item.replace(spec_orig,'')
        if 'pl(' in item:
            # Power law spectrum
            sedstart = item.find('pl(')
            sedend = item.find(')', sedstart)
            spec_orig = item[sedstart:sedend+1]
            spec = item[sedstart:sedend].replace('pl(','')
            speclist = spec.split(',')
            if speclist[0] == '1.':
                calculation['scene'][0]['spectrum']['sed'] = {'sed_type': 'powerlaw', 'index': -1*float(speclist[1]), 'unit': speclist[2]}
            item = item.replace(spec_orig,'')
        if 'bb(' in item:
            # blackbody spectrum
            sedstart = item.find('bb(')
            sedend = item.find(')', sedstart)
            spec_orig = item[sedstart:sedend+1]
            spec = item[sedstart:sedend].replace('bb(','')
            speclist = spec.split(',')
            calculation['scene'][0]['spectrum']['sed'] = {'sed_type': 'blackbody', 'teff': float(speclist[0])}
            item = item.replace(spec_orig,'')
        if 'em(' in item:
            # emission line
            sedstart = item.find('em(')
            sedend = item.find(')', sedstart)
            spec_orig = item[sedstart:sedend+1]
            spec = item[sedstart:sedend].replace('em(','')
            speclist = spec.split(',')
            line = {'center': (float(speclist[0])*u.angstrom).to_value(u.micron), 'emission_or_absorption': 'emission', 'profile': 'gaussian', \
                'strength': syn.units.convert_flux(float(speclist[0]) * u.angstrom, float(speclist[2])*syn.units.FLAM, syn.units.FNU), 'width': (float(speclist[1])/2)/float(speclist[0]) * 3e8}
            calculation['scene'][0]['spectrum']['lines'].append(line)
            item = item.replace(spec_orig,'')
        if 'rn(' in item:
            # normalize
            if 'band(' in item:
                # normalize to bandpass
                normstart = item.find('band(')
                normend = item.find(')', normstart)
                band_orig = item[normstart:normend+1]
                band = item[normstart:normend].replace('band(','')
                calculation['scene'][0]['spectrum']['normalization']['bandpass'] = band
                calculation['scene'][0]['spectrum']['normalization']['type'] = "photsys"
                item = item.replace(band_orig,'')
                del calculation['scene'][0]['spectrum']['normalization']['norm_wave']
                del calculation['scene'][0]['spectrum']['normalization']['norm_waveunit']
            if 'box(' in item:
                # normalize to wavelength
                normstart = item.find('box(')
                normend = item.find(')', normstart)
                band_orig = item[normstart:normend+1]
                band = item[normstart:normend].replace('box(','')
                bandlist = band.split(',')
                calculation['scene'][0]['spectrum']['normalization']['norm_wave'] = float(bandlist[0])
                calculation['scene'][0]['spectrum']['normalization']['norm_waveunit'] = 'angstroms'
                item = item.replace(band_orig, '')
            if 'ebmvx(' in item:
                # reddening
                redstart = item.find('ebmvx(')
                redend = item.find(')', redstart)
                red_orig = item[redstart:redend+1]
                red = item[redstart:redend].replace('ebmvx(','')
                redlist = red.split(',')
                val = float(redlist[0])
                calculation['scene'][0]['spectrum']['extinction']['bandpass'] = 'v'
                calculation['scene'][0]['spectrum']['extinction']['unit'] = 'mag'
                # TODO: We have no translations for some of these reddening laws
                if redlist[1] == "mwavg":
                    calculation['scene'][0]['spectrum']['extinction']['law'] = "mwrv31"
                    calculation['scene'][0]['spectrum']['extinction']['value'] = val*3.1
                elif redlist[1] == "mwdense":
                    calculation['scene'][0]['spectrum']['extinction']['law'] = "mwrv55" # this law is Rv = 5.5, but close enough?
                    calculation['scene'][0]['spectrum']['extinction']['value'] = val*5.0
                elif redlist[1] == "mwrv21":
                    calculation['scene'][0]['spectrum']['extinction']['law'] = ""
                    calculation['scene'][0]['spectrum']['extinction']['value'] = val*2.1
                elif redlist[1] == "mwrv4":
                    calculation['scene'][0]['spectrum']['extinction']['law'] = "mwrv40"
                    calculation['scene'][0]['spectrum']['extinction']['value'] = val*4.0
                elif redlist[1] == "lmcavg":
                    calculation['scene'][0]['spectrum']['extinction']['law'] = "lmcavg"
                    calculation['scene'][0]['spectrum']['extinction']['value'] = val*3.41
                elif redlist[1] == "lmc30dor":
                    calculation['scene'][0]['spectrum']['extinction']['law'] = ""
                    calculation['scene'][0]['spectrum']['extinction']['value'] = val*2.76
                elif redlist[1] == "smcbar":
                    calculation['scene'][0]['spectrum']['extinction']['law'] = "smcbar"
                    calculation['scene'][0]['spectrum']['extinction']['value'] = val*2.74
                elif redlist[1] == "xgalsb":
                    calculation['scene'][0]['spectrum']['extinction']['law'] = ""
                    calculation['scene'][0]['spectrum']['extinction']['value'] = val*4.05 # pulled out of the Calzetti et al. (2000. ApJ, 533, 682) paper

                calculation['calculation']['settings'] = {}
                if redend > normstart:
                    # I manually tested this, and it seems to work.  Changing the order of the reddening and band resulted in
                    #   different extinction_first values.
                    # if the reddening ended after the normalization band definition was started, we applied reddening later
                    # the normalization will have been removed when the reddening is defined, but this won't fail because the
                    # normalization entry isn't fully removed, and thus reddening couldn't start before it did.
                    calculation['calculation']['settings']['extinction_first'] = False
                else:
                    calculation['calculation']['settings']['extinction_first'] = True
                item = item.replace(red_orig,'')
            norm = item.replace(')','').split(',')
            flux = float(norm[2])
            flux_units = norm[3]
            # TODO: Revisit this unit conversion.  Pandeia supports vegamag.
            if (flux_units == 'vegamag') and (calculation['scene'][0]['spectrum']['normalization']['type'] == "at_lambda"):
                #Need to convert vegamag to mjy to normalize at_lambda
                sp = pysynphot.ArraySpectrum(np.array([calculation['scene'][0]['spectrum']['normalization']['norm_wave']]), np.array([flux]))
                sp.convert('jy')
                flux = sp.flux[0] / 1e6
                flux_units = 'mjy'
            calculation['scene'][0]['spectrum']['normalization']['norm_flux'] = flux
            calculation['scene'][0]['spectrum']['normalization']['norm_fluxunit'] = flux_units

    return calculation

def dummy_parse(adict, calculation, pweb_key):
    # dummy for now
    return calculation

def lineval_parse(adict, calculation, pweb_key):
    calculation['scene'][0]['spectrum']['lines'] = []
    if adict['target_lines'] is not None:
        for line in adict['target_lines']:
            fluxunits = FLUX_UNIT_DICT[line['fluxunits']]
            waveunits = WAVE_UNIT_DICT[line['waveunits']]

            width = ((line['fwhm'] / line['center']) * const.c.to('km/s')).value
            strength = syn.units.convert_flux(line['center'] * waveunits, line['flux']*fluxunits, syn.units.FNU)
            center = (line['center'] * waveunits).to_value(u.micron)

            newline = {'center': center, 'emission_or_absorption': 'emission', 'profile': 'gaussian', 'strength': strength, 'width': width}
            calculation['scene'][0]['spectrum']['lines'].append(newline)
    return calculation

def background_parse(adict, calculation, pweb_key):
    sky_expr = adict['sky_expr']

    # Use pysynphot to parse pyetc sky expression.
    sky_expr_parsed = SP.parse_spec(sky_expr)
    sky_expr_parsed.convert('jy')
    sky_expr_parsed.convert('micron')

    # Now assemble the spectrum (flam/arcsec and angstroms) in a format pandeia understands (MJy/sr and microns)
    bgflux = sky_expr_parsed.flux * 1e-6 / (u.arcsec * u.arcsec).to(u.sr)
    wave = sky_expr_parsed.wave

    calculation['background'] = [wave, bgflux]

    return calculation

parsefuncs = [dummy_parse, detector_parse, extraction_parse, obsmode_parse, obswave_parse, sourceval_parse, lineval_parse, background_parse]
validfuncs = [disperser_valid, simmode_valid]

H2J_INPUT_CHECKS = {
#    PENG-name,            JENG-name,                                           data-type,         can-force?,  ...
    'SNR':                 [('configuration','detector','snr'),                 'TYPE:NUM',        False, ],
    'airglow_str':         [dummy_parse,                                        'TYPE:STR',        True ],
    'aperture':            [('configuration','instrument','aperture'),          'TYPE:STR',        False ],
    'binx':                [('configuration','detector','bin_spatial'),         'TYPE:INT',        False ],
    'biny':                [('configuration','detector','bin_dispersion'),      'TYPE:INT',        False ],
    'ccddarklevel':        [('configuration','detector','dark_level'),          'TYPE:STR',        False ],
    'ccdmode':             [('strategy','ccdmode'),                             'TYPE:STR',        False ],
    'central_wavelength':  [dummy_parse,                                        'TYPE:INT',        True ],
    'crsplit':             [('configuration','detector','nsplit'),              'TYPE:INT',        False ], # see also nreads
    'detector':            [detector_parse,                                     HST_DETECTORS,     False ],
    'disperser':           [disperser_valid,                                    'TYPE:STR',        False ],
    'earthshine_str':      [dummy_parse,                                        'TYPE:STR',        True ],
    'etcid':               [dummy_parse,                                        'TYPE:STR',        False ], # will never be used here
    'extraction':          [extraction_parse,                                   'TYPE:DICT',       False ],
    'fdiameter':           [dummy_parse,                                        'TYPE:NUM',        False ],
    'filter':              [('configuration', 'instrument','filter'),           'TYPE:STR',        False ],
    'fuvglowregion':       [('configuration','detector','fuv_glow_region'),     'TYPE:STR',        True ],
    'gain':                [('configuration','detector','gain'),                'TYPE:NUM',        False ],
    'grating':             [('configuration','instrument','disperser'),         'TYPE:STR',        False ],
    'helium_str':          [dummy_parse,                                        'TYPE:STR',        True ],
    'instrument':          [('configuration','instrument','instrument'),        HST_INSTS,         False ],
    'nreads':              [('configuration','detector','nsplit'),              'TYPE:INT',        False ], # see also crsplit
    'obsmode':             [obsmode_parse,                                      'TYPE:STR',        False ],
    'mirror':              [dummy_parse,                                        HST_MIRRORS,       False ],
    'obswave':             [obswave_parse,                                      'TYPE:NUM',        False],
    'post_flash_electrons':[('configuration','detector','post_flash_electrons'),'TYPE:INT',        False ],
    'science_mode':        [('configuration','instrument','mode'),              HST_MODES,         False ],
    'sclength':            [('configuration','instrument','scan_length'),       'TYPE:NUM',        False ],
    'scrate':              [('configuration','instrument','scan_rate'),         'TYPE:NUM',        False ],
    'simmode':             [simmode_valid,                                      'TYPE:STR',        True ], # might force these for a while to pick up more stis tests
    'sky_expr':            [background_parse,                                   'TYPE:STR',        True ],
    'source_expr':         [sourceval_parse,                                    'TYPE:STR',        False ],
    'target_lines':        [lineval_parse,                                      'TYPE:LIST',       True ],
    'time':                [('configuration','detector','time'),                'TYPE:NUM',        False ],
    'zodi_str':            [dummy_parse,                                        'TYPE:STR',        True ],
}

def pyetc_in_dict_to_pandeia_in_dict(adict, input_fname, pweb_dict=None):
    # init
    new_inputs = {}
    msgs = []

    # first pass - go through and warn about unused args/values
    for key in adict:
        # basics
        assert key in H2J_INPUT_CHECKS, 'UNSUPPORTED (as yet!) pyetc input arg: %s' % key
        arg_info = H2J_INPUT_CHECKS[key]
        value = adict[key]
        newname = arg_info[0] if arg_info[0] else key
        valid_vals = arg_info[1]
        assert not valid_vals is None, 'Error in data defined in H2J_INPUT_CHECKS for: %s' % key
        if valid_vals == 'TYPE:DICT':
            assert isinstance(value, dict) or value is None, 'Expected type dict for value of: %s' % key
            # no more checks here yet
            new_inputs[newname] = value
        elif valid_vals == 'TYPE:LIST':
            assert isinstance(value, list) or value is None, 'Expected type list for value of: %s' % key
            # no more checks here yet
            new_inputs[newname] = value
        elif valid_vals == 'TYPE:STR':
            assert isinstance(value, str) or value is None, 'Expected type str for value of: %s' % key
            # no more checks here yet
            if value is not None:
                new_inputs[newname] = value.replace('\n','')
            else:
                new_inputs[newname] = value
        elif valid_vals == 'TYPE:INT':
            assert type(value) == int, 'Expected type int for value of: %s' % key
            # no more checks here yet
            new_inputs[newname] = value
        elif valid_vals == 'TYPE:NUM':
            assert type(value) in (int, float), 'Expected type int or float for value of: %s' % key
            # no more checks here yet
            new_inputs[newname] = value
        else:
            assert type(valid_vals) == list and len(valid_vals) > 0, 'Error in data defined in H2J_INPUT_CHECKS for: %s' % key
            assert value in valid_vals, f'Error in defined data: {newname} {value}'
            new_inputs[newname] = value

    # instrument and mode, filters and dispersers
    calculation = build_default_calc("jwst", "nircam", "sw_imaging")
    calculation['background'] = {}
    # some calculations have no filter or disperser
    calculation['configuration']['instrument']['filter'] = None
    calculation['configuration']['instrument']['disperser'] = None
    calculation['configuration']['instrument']['aperture'] = None
    calculation['configuration']['detector'] = {} # Detector is full of JWST-only keys we'll replace.
    del calculation['background_level']
    for item in new_inputs:
        if isinstance(item, tuple):
            if len(item) == 2:
                calculation[item[0]].update({item[1]: new_inputs[item]})
            elif len(item) == 3:
                calculation[item[0]][item[1]].update({item[2]: new_inputs[item]})
        elif isinstance(item, str):
            calculation.update({item: new_inputs[item]})
    for item in new_inputs:
        if item in parsefuncs:
            if pweb_dict:
                calculation = item(adict, calculation, pweb_dict)
            else:
                calculation = item(adict, calculation, input_fname)
    for item in new_inputs:
        if item in validfuncs:
            calculation = item(adict, calculation)

    # merge both pyetc wfc3 instruments to wfc3
    if "wfc3" in calculation['configuration']['instrument']['instrument']:
        calculation['configuration']['instrument']['instrument'] = "wfc3"

    return calculation

def write_sorted_json(data, fname):
    with open(fname, 'w') as f:
        json.dump(data, f, indent=4, separators=(',', ': '), sort_keys=True, cls=NumPyArangeEncoder)
        f.write('\n')


def convert_dict(peng_dict, input_fname):
    pandeia_input_dict = pyetc_in_dict_to_pandeia_in_dict(peng_dict, input_fname)
    return pandeia_input_dict


def convert_file(input_fname, top_out_dir=None, limit_outputs=None):
    """ Convert a single PENG file to a JENG file """
    # output
    # names and paths
    path_part, fname_part = os.path.split(input_fname)
    assert path_part.startswith('/'), 'unexpected relativity in input fname: %s' % input_fname # not absolute
    basename = os.path.splitext(fname_part)[0]
    outstr = '\n'+basename+' ...\n'
    outstr += '-'*80; outstr += '\n'
    out_fname = basename+'.jeng' # for dumping it locally
    # top dir specified for the ouput?
    if top_out_dir:
        unused, path_part_btm = path_part.split('/pyetc/test/')
        path_part = top_out_dir+os.sep+path_part_btm
    if not os.path.exists(path_part):
        os.makedirs(path_part)

    # output name
    out_fname = path_part+os.sep+basename+'.jeng' # for dumping it in place
    # load it to a dict
    pyetc_input_dict = peng2dict(input_fname)
    # convert it
    pandeia_input_dict = convert_dict(pyetc_input_dict, input_fname)


    instdotmode = '%s.%s' % (pandeia_input_dict['configuration']['instrument']['instrument'], pandeia_input_dict['configuration']['instrument']['mode'])

    # see if they asked us to limit the number of test files we create
    if limit_outputs:
        assert limit_outputs.get(instdotmode, 0) < 500, 'Have enough %s tests now' % instdotmode

    # write it out
    if os.path.exists(out_fname):
        os.remove(out_fname)
    write_sorted_json(pandeia_input_dict, out_fname)
    outstr += 'Wrote %s test file: %s\n' % (instdotmode, out_fname)
    return outstr


def run_some_samples(pyetc_test_dir, top_out_dir, single_file=None):
    """
    Select a sampling of all the modes we can run.
    Either send in the first 2 args, or just use the single_file arg.
    This could be smoother, I know.
    """

    # paths under pyetc test
    # places = [
    #     'engine/test_kits/geom-detector', # put these in first so they get grabbed
    #     'engine/bit/stis/imag', # about 250 tests
    #     'engine/bit/wfc3ir/imag/ir',
    #     'engine/sample', # has some COS (not done yet) but some STIS that might work
    #     'engine/spider/acs/imag/sbc-stars-1',
    #     'engine/spider/common/wfc3ir/imag', # this one line adds about 1300 tests
    #     'engine/spider/stis/spec/stellar-pt',
    #     'engine/spider/stis/spec/stellar-pt-bin',
    #     'engine/spider/wfc3ir/imag', # this includes both pt and ext srcs
    #     'engine/spider/wfc3uvis/imag', # this includes both pt and ext srcs
    #     'engine/spider/wfc3ir/scimag',
    #     'engine/test_kits/geom-sim',
    #     'engine/test_kits/limits',
    # ]
    places = [
        'engine/test_kits/', # put these in first so they get grabbed
        'engine/bit/',
        'engine/sample', # has some COS (not done yet) but some STIS that might work
        'engine/spider/'
    ]
    flist = []
    skips = {}

    if single_file:
        flist = [single_file]
        top_out_dir = None # will drop new jeng file right next to old peng file
    else:
        for place in places:
            sublist = utils.rglob(pyetc_test_dir+'/'+place, '*.peng')
            assert len(sublist) > 0, 'no pengs found under %s ??' % place
            flist += sublist

    done = {}
    for fname in flist:
        # outstr = convert_file(fname, top_out_dir)#, limit_outputs=done)
        # # if no exception than print the ouput
        # print(outstr)
        # instdotmode = get_instdotmode_from_output(outstr)
        # if instdotmode in done:
        #     done[instdotmode] += 1
        # else:
        #     done[instdotmode] = 1
        try:
            outstr = convert_file(fname, top_out_dir)#, limit_outputs=done)
            # if no exception than print the ouput
            instdotmode = get_instdotmode_from_output(outstr)
            if instdotmode in done:
                done[instdotmode] += 1
            else:
                done[instdotmode] = 1
        except Exception as e:
            if len(flist) == 1:
                raise e
            reason = str(e) if type(e) == AssertionError else repr(e)
            if reason in skips:
                skips[reason].append(fname)
            else:
                skips[reason] = [fname]

    total_skipped = 0
    if skips:
        print("\n\nSKIPS:\n")
    for reason in skips:
        ncases = len(skips[reason])
        print("Skipped %d cases of: %s,  e.g.: %s" % (ncases, reason, skips[reason][0]))
        total_skipped += ncases

    print('\nConverted %d and skipped %d tests. Of the converted: %s' % (len(flist)-total_skipped, total_skipped, done))


#
# main routine
#
if __name__=='__main__': # in case something else imports this file
    # If this is run with no args, go fetch a same set of tests
    if len(sys.argv) == 1:
        run_some_samples(os.path.expanduser(f'{DATA_DIR}/../pyetc/test'),
                         os.path.expanduser(f'{DATA_DIR}/../pandeia_test/tests/engine/hst'))

    # If this is given a peng file name as an arg on the command line, translate it
    else:
        for fname in sys.argv[1:]:
            run_some_samples(None, None, single_file=fname)

    sys.exit(0)
