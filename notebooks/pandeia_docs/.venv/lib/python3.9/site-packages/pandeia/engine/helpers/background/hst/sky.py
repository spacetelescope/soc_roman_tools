from __future__ import division

import datetime
import os

from pandeia.engine.helpers.background.hst import util
from pandeia.engine.helpers.background.hst import sky_angles as angles
from pandeia.engine.custom_exceptions import BackgroundError, RangeError, UserInputError

from pandeia.engine.utils import recursive_subclasses, merge_data

'''
keep the synphot expression generation section of the code in sync with 
https://github.com/spacetelescope/pyetc/blob/master/pyetc/etc_web/etc/sky.py
'''
#TODO: sky.dat needs to migrate to CDBS & have a versioned name
MODELS = util.read_dict(os.path.join(os.path.dirname(__file__),
                                'sky.dat'))

class SkyElement(object):
    """Base class for actual elements of the sky model"""
    #Subclasses will override this value
    formkeys = []


    def __init__(self, formdata=None):
        """formdata=a dictionary containing the key/value pairs required to
        instantiate this SkyElement
        """
        self.fd=dict()
        self.descr='Not yet set'
        self.expr='Not yet set'

        if formdata:
            #pick out what we care about
            for k in self.formkeys:
                self.fd[k]=formdata[k]
                #This will raise a KeyError if any keys of interest
                #are not present. The caller (typically the Sky constructor)
                #will have to deal with this exception

    def __str__(self):
        return str(self.descr)


    def make_synexpr(self):
        """Subclasses must implement this method. It should set both the
        .expr and the .descr attributes
        """
        raise NotImplementedError("Sorry, this must be implemented in the subclass")


class Airglow(SkyElement):
    """Airglow is modeled as a set of narrow emission lines in the UV."""

    #The only input values we care about
    formkeys = ['standard']
    lines = ['el1215a',
             'el1302a',
             'el1356a',
             'el2471a']

    #Normalization constants for each level
    AIRGLOW_LEVELS = {
        'none': None,
        'low': {'el1215a': 0.20,
                'el1302a': 0.01333333333,
                'el1356a': 0.012,
                'el2471a': 0.010},
        'average': 1.0,
        'high': 2.0
        }


    def make_synexpr(self):
        self.descr=self.fd['standard']
        if self.descr == 'none':
            self.expr = None
            return

        elif self.descr == 'low':
            elements=[]
            #Special case, different multipliers for each line
            norm = self.AIRGLOW_LEVELS[self.descr]
            for l in self.lines:
                elements.append("%f*spec(%s)"%(norm[l],MODELS[l]))
            self.expr = ' + '.join(elements)

        else: #All others have a single multiplier applied to all lines
            elements = []
            for l in self.lines:
                elements.append("spec(%s)"%MODELS[l])
            self.expr = '(%s)'%' + '.join(elements)
            if self.AIRGLOW_LEVELS[self.descr] != 1.0:
                self.expr = "%f*(%s)"%(self.AIRGLOW_LEVELS[self.descr], self.expr)


class Earthshine(SkyElement):
    """Earthshine is modelled by a spectrum file."""

    formkeys = ['standard','magnitude','multiplier','choice']
    specfile = "spec(%s)"%MODELS['earthshine']
    # These are scalar multipliers
    EARTHSHINE_LEVELS = {
        'none': None,
        'low': None,
        'average': 0.5,
        'high': 1.0,
        'extremely_high': 2.
        }

    def __init__(self, formdata=None):
        #Call the parent
        SkyElement.__init__(self, formdata)
        #This tells us which type of calculation we're going to do
        self.etype = self.fd['choice']
        self._set_norm()


    def make_synexpr(self):
        if self.norm is None:
            self.expr = None
        else:
            if self.etype in ('magnitude'):
                self.expr = "rn(%s,band(johnson,v),%f,vegamag)"%(self.specfile,self.norm)
            elif self.etype in ('standard', 'multiplier'):
                self.expr = self.specfile
                if self.norm != 1.0:
                    self.expr = "%f*%s"%(self.norm, self.expr)
            else:
                raise UserInputError(f'Unrecognized normalization option {self.etype}')


    def _set_norm(self):
        if self.etype == 'standard':
            level = self.fd[self.etype]
            self.norm = self.EARTHSHINE_LEVELS[level]
            self.descr = level

        elif self.etype == 'magnitude':
            self.norm = self.fd[self.etype]
            self.descr = 'User Defined (normalized to %f V)'%self.norm

        elif self.etype == 'multiplier':
            self.norm = self.fd[self.etype]
            self.descr = 'User Defined (scaled by %f)'%self.norm

        else:
            raise UserInputError(f'Unrecognized normalization option {self.etype}')


class HeliumEmiss(SkyElement):
    """HeliumEmiss is modelled by a spectrum file."""

    formkeys = ['standard']
    specfile = "spec(%s)"%MODELS['helium']

    # These are scalar multipliers
    HELIUM_EMISSION_LEVELS = {
        'none':      None,
        'average':   0.1,
        'high':      0.5,
        'very high': 1.5,
        }

    def __init__(self, formdata=None):
        SkyElement.__init__(self, formdata)
        self.descr = self.fd['standard']
        self.norm = self.HELIUM_EMISSION_LEVELS[self.descr]

    def make_synexpr(self):
        if self.norm is None:
            self.expr = None
        else:
            self.expr = self.specfile
            if self.norm != 1.0:
                self.expr = "%f*%s"%(self.norm, self.expr)


class ZodiacalLight(SkyElement):
    """Zodiacal light is modeled by a spectrum file that can be
    normalized in a variety of ways, including based on the angle
    between the sun and the observed coordinates.
    """
    #TODO - Although only a subset of these keys are necessary to
    #construct a valid ZodiacalLight object, at present all are
    #required. We might override the __init__ completely and make it
    #complicated enough to cope with conditionals here.
    formkeys = ['choice','standard','ra','dec',
                'sunattribute', 'sun', 'magnitude', 'multiplier',
                'date', 'month', 'day']
    specfile = "spec(%s)"%MODELS['zodi']
    ZODI_LEVELS = {  #Normalized to Johnson V
        'none': None,
        'low': 23.3,
        'average': 22.7,
        'high': 22.1
    }

    #TODO: handle the posfile in the same way as the other data files
    #  except that some variable expansion has to be done, unlike for
    #  the other files for which that is handled by pysynphot
    #posfile = MODELS['zodipos']

    def __init__(self, formdata=None):
        #Then call the parent
        SkyElement.__init__(self, formdata)
        #This tells us which type of calculation we're going to do
        self.ztype = self.fd['choice']
        self._set_norm()

    def make_synexpr(self):
        if self.norm is None:
            self.expr = None
        else:
            if self.ztype in ('standard', 'magnitude', 'coords'):
                self.expr = "rn(%s,band(johnson,v),%f,vegamag)"%(self.specfile,self.norm)
            elif self.ztype in ('multiplier'):
                self.expr = "%f*%s"%(self.norm,self.specfile)

    def _set_norm(self):
        #There are many ways to specify how the Zodiacal Light spectrum should be normalized
        if self.ztype == 'standard':
            level = self.fd[self.ztype]
            self.norm = self.ZODI_LEVELS[level]
            self.descr = level

        elif self.ztype == 'magnitude':
            self.norm = self.fd[self.ztype]
            self.descr = 'User Defined (normalized to %f V)'%self.norm

        elif self.ztype == 'multiplier':
            self.norm = self.fd[self.ztype]
            self.descr = 'User Defined (scaled by %f)'%self.norm

        elif self.ztype == 'coords':
            zra, zdec = self.fd['ra'], self.fd['dec']
            ra, dec = angles.radec_to_decimal( zra, zdec )

            if self.fd['sunattribute'] == 'date':
                year, month, day = self.fd['date'], self.fd['month'], self.fd['day']
                date = datetime.date(year,month,day)
                helong, helat, sunangle = angles.radec_date_to_helonglat(ra, dec, date)
                self.descr = ('%s %s %s/%s/%s (%f deg.)' % (zra,zdec,month,day,year,sunangle))

            elif self.fd['sunattribute'] == 'sunangle':
                zsun = self.fd['sun']
                # Table is defined from 0-180, but users can put in anything
                sunangle = angles.limit_longitude(float(zsun))
                self.descr = ('%s %s Sun Angle = %f deg.' % (zra, zdec, zsun))
                helong, helat = angles.radec_sunangle_to_helonglat(ra, dec, sunangle)

            else:
                raise UserInputError('Unknown sun specification %s'%self.fd['sunattribute'])
            try:
                self.norm = angles.zodi_lookup(helong, helat)
            except RangeError as e:
                e.addinfo("These values are derived from the RA, Dec, and Date or Helio-ecliptic Longitude. You must either select different values for these quantities, or a different normalization method for the zodiacal light.")
                raise

        else:
            raise UserInputError(f'Unknown normalization type {self.ztype}')


class Sky(object):
    """Class to compute the sky background synphot expression including all relevant SkyElements.
    Presently supports these elements:
        Airglow, Earthshine, and ZodiacalLight.
    Not all the elements are required. Elements will be picked out based on
    the contents of the dictionary that is passed in."""
    default_config = {
        'zodiacallight': {'choice': 'standard', # standard, coords, magnitude, multiplier
                 'standard': 'average',
                 'ra': '00:00:00',
                 'dec': '00:00:00',
                 'sunattribute': 'date', # date, sunangle
                 'date': 2009,
                 'month': 1,
                 'day': 1,
                 'sun': 90.0,
                 'magnitude': 30.0,
                 'multiplier': 1.0},
        'earthshine': {'choice': 'standard',  # standard, magnitude, multiplier
                       'standard': 'average',
                       'multiplier': 1.0, 
                       'magnitude': 30.0},
        'airglow': {'choice': 'standard',
                    'standard': 'average'} # standard
        }

    def __init__(self, config={}, webapp=False, **kwargs):
        self.webapp = webapp

        self.expr = None
        types = recursive_subclasses(SkyElement)
        elements = [t.__name__.lower() for t in types]
        self.sky_map = dict(list(zip(elements, types)))
        self.bg_input = merge_data(self.default_config, config)

    def make_synexpr(self):
        elements = []
        for key in sorted(self.bg_input.keys()):
            sky_elem = self.sky_map[key](self.bg_input[key])
            sky_elem.make_synexpr()
            if sky_elem.expr is not None:
                elements.append(sky_elem.expr)

        if len(elements) > 0:
            self.expr = ' + '.join(elements)

