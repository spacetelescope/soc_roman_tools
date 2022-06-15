'''Map the source form parameters into a form that the etc_engine can use.

Two approaches are used. The old is to construct a synphot source
expression. The new is to construct a pysynphot spectrum object.
Both take a form dictionary as input.
'''
from __future__ import division
import os, sys

import pandeia.engine.helpers.bit.pyetc_util as util

#Set up some constants
wave_scaling = {'microns':    10.0**4,
                'nanometers': 10.0,
                'angstroms':   1.0}
symbol = {'microns': "&mu;",
          'nanometers': "nm",
          'angstroms': "&Aring"}

# Map of target choices to synphot files/expressions
SDICT = util.read_dict(os.path.join(os.path.dirname(__file__),'etc_web_data/targets.dat'))

# Map target kind onto target value selector variable
SPECNAMEMAP = {
    'SpectrumCK':'fSpectrumCK',
    'SpectrumPickles':'fSpectrumPickles',
    'SpectrumKurucz':'fbpgsfile',
    'SpectrumStellar':'fStellar',
    'SpectrumNonStellar':'fnonstellar',
    'SpectrumHST': 'fcalfile',
    'SpectrumPhoenix': 'fphoenixfile'
    }

# TARGET_HELPERS massage the raw data of the target mapping file into description and
# pysynphot expression tuples,  basically user results text + psyn code.
TARGET_HELPERS = {
    # form_var    :  lambda key, target.dat[form_var][key] :   description,  psyn_expr
    "fSpectrumCK"      : lambda k, v:  ("Castelli-Kurucz Models " + k + " " + v,      'icat(ck04models,%s,0.0,%s)' % tuple(v.split())),
    "fbpgsfile"        : lambda k, v:  ("Kurucz Models " + k + " " + v,               'icat(k93models,%s,0.0,%s)' % tuple(v.split())),
    "fSpectrumPickles" : lambda k, v:  ("Pickles Models " + v[0],                         v[1]),
    "fcalfile"         : lambda k, v:  ("Standard HST Star Spectra " + v[0],          "spec(%s)" % (v[1],)),
    "fphoenixfile"     : lambda k, v:  ("Phoenix M dwarf Models " + v[0],             "spec(%s)" % (v[1],)),
    "fStellar"         : lambda k, v:  ("Bruzual Synthetic Stellar Spectra " + v[0],  "spec(%s)" % v[1]),
    "fnonstellar"      : lambda k, v:  ("Non-stellar Objects " + v[0],                "spec(%s)" % v[1]),
}

# Provides filter names for use in building output page. Neither .capitalize()
# or .upper() work as a generic solution; we need to apply one of them and
# counteract exception cases with this dictionary.
filter_output_names = {"SDSS":   "Sloan",
                       "ACS":    "ACS",
                       "NICMOS": "NICMOS",
                       "WFC3":   "WFC3"}

def map_filter(filter):
    '''convert the web form values for the normalization bandpass into a
    form that synphot/pysynphot can use'''
    #Normally, we can just replace a / by a comma and we're done.
    filter = filter.lower().replace('/',',')

    #But we have to change sloan to sdss
    if filter.startswith('sloan'):
        filter = filter.replace('sloan','sdss')

    #and specify the detector for NICMOS.
    #JETC uses nicmos, 2: ETC/datafiles/SynphotBands.xml, lines 110:123
    elif filter.startswith('nicmos'):
        filter = filter.replace('nicmos','nicmos,2')

    # and for WFC3UVIS
    elif filter.startswith('wfc3'):
        filter = filter.replace('uvis','uvis1')

    # the ACS detector is not specified; we arbitrarily pick WFC
    elif filter.startswith('acs'):
        filter = filter.replace('acs','acs,wfc1')
    return filter

def map_name(stype, fname):
    '''determine specific file matching given name or key

    >>> map_name('SpectrumPickles', 'O5V')
    ('Pickles O5V', 'crgrid$pickles/dat_uvk/pickles_uk_1.fits')
    '''
    try:                              # .e.g. stype=SpectrumPickles
        keyvar = SPECNAMEMAP[stype]   # e.g.  keyvar=fSpectrumPickles
        return TARGET_HELPERS[keyvar](fname, SDICT[keyvar][fname])
    except KeyError:
        raise ValueError("Unrecognized source type: %s, %s"%(stype,fname))

class Source(object):
    """
    This class embodies the source descriptions and abiltiy to construct
    pysynphot spectrum objects.
    """
    def __init__(self, request_dir, formdata=None, obsmode=None):
        """Construct source object"""
        self.request_dir = request_dir
        self.obsmode=obsmode
        if formdata:
            self.inputs_from_form(formdata)
        else:
            self.expr = ''
            self.help = '' #for output web form

    def inputs_from_form(self, formdata):
        '''Construct a synphot source expression from the form info.
        Currently forms consist of mostly exclusive options in dict
        format'''

        # These functions construct both the synphto expression element
        # and the user-friendly string displayed on web output.

        # specexpr is the fundamental source expression, which can
        # be redshifted, extincted, renormalized, and had emission lines
        # added to it. The order of these operations matters.
        spechelp, specexpr = self.handle_source(formdata)
        ehelp, estr        = self.handle_extinction(formdata)
        rnhelp, rnstr      = self.handle_renormalization(formdata,obsmode=self.obsmode)
        zhelp, zstr     = self.handle_redshift(formdata)
        lhelp, lines    = self.handle_emission_lines(formdata)

        #First apply redshift
        specexpr = zstr % specexpr

        #Then apply extinction either before or after normalization
        if formdata['fextinctiontype'] == 'before':
            specexpr += estr
            specexpr = rnstr % specexpr
        else:
            specexpr = rnstr % specexpr
            specexpr += estr

        #Finally, assign to the expr attribute.
        #Lines will be passed in separately.
        self.expr = specexpr



        #Also construct the "help" attribute for user friendliness.
        #Note emission lines returns non-standard multi-line help

        #TODO: Make pysyn_help display smarter so it only shows when
        #mode=debug; or, filter out filenames. In the meantime, disable
        #this display for security reasons.
        pysyn_help = ["Pysynphot:", self.expr]
        self.help = [spechelp, ehelp, rnhelp, zhelp] + lhelp #+ [["<p/>",""], pysyn_help]

        #Verify that we got a result out of all that work before returning.
        if not self.expr:
            #Lines might be specified separately
            if lines:
                #Keep pysynphot parser from failing: specify flat spectrum
                #that is uniformly zero
                self.expr = 'unit(0.0,flam)'
            else:
                raise ValueError('Source form selections did not result in a valid source specification.\n These might include no flux from the selected source combination.')

    def handle_source(self, fd):
        # Always uses key: fsorigin
        # Conditionally uses one of:
        #  fSpectrumCK, fSpectrumPickles, fbpgsfile, fStellar, fcalfile, fphoenixfile, fnonstellar,fbbtemp,(fIndex,fIsLambda)
        # determine which option to use
        stype = fd['fsorigin']
        if stype == 'SpectrumUpload':
            specexpr = fd['fUploadFileLocal']
            if self.request_dir not in [".","./"]:
                specexpr = self.request_dir + '/' + specexpr
            info = "User uploaded spectrum:  " + fd['fUploadFile']
        elif stype == 'SpectrumHstOther':
            specexpr = fd['fOtherFile']
            info = "Other HST spectrum " + repr(fd['fOtherFile'])
        elif stype in ['SpectrumCK', 'SpectrumPickles','SpectrumKurucz','SpectrumStellar',
                       'SpectrumNonStellar','SpectrumHST', 'SpectrumPhoenix']:
            selected = fd[SPECNAMEMAP[stype]]
            info, specexpr = map_name(stype, selected)

        #The following options have sanity checking applied to them:
        elif stype == 'SpectrumBlackBody':
            kelvins = float(fd['fbbtemp'])
            if kelvins <= 0:
                raise ValueError("Temperature %f must be positive"%kelvins)
            specexpr = 'bb(%f)' % kelvins
            info = "Black body at " + str(kelvins) + "K"

        elif stype == 'SpectrumPowerLaw':
            specexpr = 'pl(1.,%f,flam)' % fd['fIndex']
            info = "Power Law with index " + str(fd['fIndex'])
        elif stype == 'SpectrumFlat':
            units = ['flam','fnu'][int(fd['fIsLambda']=='false')]
            specexpr = 'unit(1.,%s)' % units
            info = "Flat in units of " + {"flam" : "F&lambda;",
                                          "fnu"  : "F&nu;", } [units]
        elif stype == 'SpectrumEmpty':
            specexpr = ''
            info = "Empty"
        else:
            raise ValueError({},'Unknown option %s for source spectrum'%stype)
        return ["Spectrum:",info], specexpr

    EXTINCTION_MAP  = {
        "mwavg"   : "Milky Way Diffuse (Rv=3.1)",
        "mwdense" : "Milky Way Dense (Rv=5.0)",
        "mwrv21"  : "Milky Way CCM (Rv=2.1)",
        "mwrv4"   : "Milky Way CCM (Rv=4.0)",
        "lmcavg"  : "LMC Diffuse (Rv=3.41)",
        "lmc30dor": "LMC 30DorShell (Rv=2.76)",
        "smcbar"  : "SMC Bar (Rv=2.74)",
        "xgalsb"  : "Starburst (attenuation law)",
        }

    def handle_extinction(self, fd):
        # extinction model section
        #Always uses key: febv
        #Conditionally uses key: febmvtype
        ebv = fd['febv']
        if ebv:
            etype = fd['febmvtype']
            estr = ' * ebmvx(%f,%s)' % (ebv, etype)
            ehelp = self.EXTINCTION_MAP[ etype ] + " = " + str(ebv)
            ehelp += " applied " + fd['fextinctiontype'] + " normalization"
        else:
            estr = ''
            ehelp = 'None'
        return ["Extinction E(B-V):", ehelp], estr

    def handle_redshift(self, fd):
        # redshift handling
        # Always uses key: fRedshift
        z = fd['fRedshift']
        if float(z) < 0:
            raise ValueError("Redshift %f must be positive"%z)
        if z:
            zstr = 'z(%s' + ',%f)' % z
            zhelp = "z=" + str(z)
        else:
            zstr = '%s'
            zhelp = "None"
        return ["Redshift:", zhelp], zstr

    def _build_bandpass(self, fd):

        instrument = fd['instrument']
        science_mode = fd['science_mode']

        module_name = 'pandeia.engine.helpers.bit.instruments.hst.%s.web.%s_web' % (instrument, instrument)
        make_obsmode_function_name = 'make_%s_%s_obsmode' % (instrument, science_mode)

        module = util.dynamic_import(module_name)
        make_obsmode_function = getattr(module, make_obsmode_function_name)

        # STIS needs special handling: the interface to the function
        # that builds the obsmode string is non-standard and depends
        # on the science mode.
        if instrument == 'stis':
            # we delegate to the instrument_from_form() method the
            # details on how to process the input in a way that
            # makes sense to the obsmode building function.
            instrument_from_form_function_name = 'instrument_from_form'
            instrument_from_form_function = getattr(module, instrument_from_form_function_name)
            idict  = instrument_from_form_function(science_mode, fd)

            if science_mode == 'imaging':
                results = make_obsmode_function(idict['detector'], idict['aperture'])
            elif science_mode == 'targetacquisition':
                results = make_obsmode_function(idict['detector'],
                                                idict['aperture'],
                                                idict['ccdmode'],
                                                (fd['fsourceType'] == 'point'))
            else: # spectroscopic
                results = make_obsmode_function(fd)

        else:
            results = make_obsmode_function(fd)

        # the XXX_web.py API has a bit of internal variation....
        if type(results) == str:
            return results
        else:
            return results[0]

    def handle_renormalization(self, fd, obsmode=None):
        # renormalization section
        # Always uses key: fnormalize
        # Conditionally uses key: fftype_filters, rn_flux_bandpass, rn_flux_bandpass_units,
        rntype = fd['fftype']
        # What we're doing here is building a string into which the
        # spec_expr will be substituted. This is to support the option
        # of applying the extinction either before or after the
        # renormalization.
        stype = fd['fsorigin']

        if (rntype == 'fno' or stype == 'SpectrumEmpty'): #cannot renormalize nothing
            rnstr = '%s'
            rnhelp = "None"

        elif rntype.startswith('fnormalize'):

            # pick off the normalization mode from the rntype string
            junk, norm_mode = rntype.split('.',1)

            # normalize at a specific wavelength.
            if norm_mode == 'lambda':
                fflux = fd['rn_flux_lambda']
                filter = fd['rn_flux_lambda_units']

                flambda = fd['rn_lambda']
                lambda_units = fd['rn_lambda_units']
                flambda_angstrom = flambda * wave_scaling[lambda_units]

                # Bandpass is box 1.0 A wide.
                rnstr = 'rn(%s,' +'box(%f,1.),%g,%s)' % (flambda_angstrom, fflux, filter)
                rnhelp = "Flux " + str(fflux) + " " + str(filter) + " at " + \
                         str(flambda) + " " + symbol[lambda_units]

            # normalize by filter or bandpass
            else:
                fflux = fd['rn_flux_bandpass']
                funits = fd['rn_flux_bandpass_units']
                filter_type = fd['fftype_filters'].split('.',1)[1]

                if filter_type == 'observation':
                    # normalize by bandpass used for observation
                    fband = self._build_bandpass(fd)
                    band_type = " bandpass defined above."

                else:
                    # normalize by user-chosen filter
                    filter_selector = "filter." + filter_type
                    fband = map_filter(fd[filter_selector])

                    # good-looking output is tricky.
                    fband_output = fband.upper().replace(',', '/')
                    tokens = fband_output.split('/')
                    tokens[0] = filter_output_names.get(tokens[0], tokens[0].capitalize())
                    band_type = " filter " + '/'.join(tokens)

                rnstr = 'rn(%s,' + 'band(%s),%g,%s)' % (fband, fflux, funits)
                rnhelp = "Renormalized to " + funits + " = " + str(fflux) + " in " + band_type

        else:
            raise ValueError('Unrecognized normalize option %s'%rntype)
        return ["Normalization:", rnhelp], rnstr

    def handle_emission_lines(self, fd):
        """ Produce a list of dicts that specify the lines.
        See engine/observable.ObservableItem docstring for
        details.
        """
        # emission line section

        #Central wavelength may be either in microns or angstroms.
        #FWHM is always in angstroms
        #We will convert both to angstroms for pysynphot.
        self.lines = list()
        waveunits = fd.get('wavelengthUnits','angstroms')
        wavesymbol = symbol[waveunits]
        scl = wave_scaling[waveunits]
        info = [["Emission Lines:",""]]
        #Up to three lines
        for n in (1,2,3):
            c, w, f = fd['fL%d_center'%n], fd['fL%d_fwhm'%n], fd['fL%d_flux'%n]
            if c*w*f > 0:
                self.lines.append(dict(center = c*scl,
                                       fwhm = w,
                                       flux = f,
                                       waveunits = 'angstroms',
                                       fluxunits = 'flam')
                                  )
                info.append(["Line at", "Center " + str(c) + " " + wavesymbol + "&nbsp;&nbsp;&nbsp; "
                             "FWHM " + str(w) + " &Aring;&nbsp;&nbsp;&nbsp;"
                             "Flux " + str(f) + " F&lambda;"])

        if not self.lines:
            self.lines = None
            info = [["Emission Lines:", "None"]]

        return [info, self.lines]
