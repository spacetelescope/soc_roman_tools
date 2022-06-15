from __future__ import division

import datetime, os

from pandeia.engine.helpers.background.hst.util import read_dict
from pandeia.engine.helpers.background.hst import sky_angles as angles
from pandeia.engine.custom_exceptions import RangeError

'''
keep the synphot expression generation section of the code in sync with 
https://github.com/spacetelescope/pyetc/blob/master/pyetc/etc_web/etc/sky.py

keep everything except Sky class, formkeys, fd, descr, etype, ztype, level keys and error types in sync with 
https://github.com/spacetelescope/pandeia/blob/master/engine/pandeia/engine/helpers/background/hst/sky.py
'''
#TODO: sky.dat needs to migrate to CDBS & have a versioned name
MODELS = read_dict(os.path.join(os.path.dirname(__file__), '../background/hst/sky.dat'))

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
    formkeys = ['AirglowStandard']
    lines = ['el1215a',
             'el1302a',
             'el1356a',
             'el2471a']

    #Normalization constants for each level
    AIRGLOW_LEVELS = {
        'None': None,
        'Low': {'el1215a': 0.20,
                'el1302a': 0.01333333333,
                'el1356a': 0.012,
                'el2471a': 0.010},
        'Average': 1.0,
        'High': 2.0
        }


    def make_synexpr(self):
        self.descr=self.fd['AirglowStandard']
        if self.descr == 'None':
            self.expr = None
            return

        elif self.descr == 'Low':
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

    formkeys = ['EarthshineStandard','EarthshineMag','EarthshineMult','EarthshineSpec']
    specfile = "spec(%s)"%MODELS['earthshine']
    # These are scalar multipliers
    EARTHSHINE_LEVELS = {
        'None': None,
        'Low': None,
        'Average': 0.5,
        'High': 1.0,
        'Extremely_high': 2.
        }

    def __init__(self, formdata=None):
        #Call the parent
        SkyElement.__init__(self, formdata)
        #This tells us which type of calculation we're going to do
        self.etype = self.fd['EarthshineSpec']
        self._set_norm()


    def make_synexpr(self):
        if self.norm is None:
            self.expr = None
        else:
            if self.etype in ('EarthshineMag'):
                self.expr = "rn(%s,band(johnson,v),%f,vegamag)"%(self.specfile,self.norm)
            elif self.etype in ('EarthshineStandard', 'EarthshineMult'):
                self.expr = self.specfile
                if self.norm != 1.0:
                    self.expr = "%f*%s"%(self.norm, self.expr)
            else:
                raise ValueError("Unrecognized normalization option %s"%self.etype)

    def _set_norm(self):
        if self.etype == 'EarthshineStandard':
            level = self.fd[self.etype]
            self.norm = self.EARTHSHINE_LEVELS[level]
            self.descr = level

        elif self.etype == 'EarthshineMag':
            self.norm = self.fd[self.etype]
            self.descr = 'User Defined (normalized to %f V)'%self.norm

        elif self.etype == 'EarthshineMult':
            self.norm = self.fd[self.etype]
            self.descr = 'User Defined (scaled by %f)'%self.norm

        else:
            raise ValueError("Unrecognized normalization option %s"%self.etype)


class HeliumEmiss(SkyElement):
    """HeliumEmiss is modelled by a spectrum file."""

    formkeys = ['HeliumStandard']
    specfile = "spec(%s)"%MODELS['helium']

    # These are scalar multipliers
    HELIUM_EMISSION_LEVELS = {
        'None':      None,
        'Average':   0.1,
        'High':      0.5,
        'Very High': 1.5,
        }

    def __init__(self, formdata=None):
        SkyElement.__init__(self, formdata)
        self.descr = self.fd['HeliumStandard']
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
    formkeys = ['ZodiSpec','ZodiStandard','ZodiRA','ZodiDec',
                'ZodiSunAttribute', 'ZodiSun', 'ZodiMag', 'ZodiMult',
                'ZodiDate', 'ZodiMonth', 'ZodiDay']
    specfile = "spec(%s)"%MODELS['zodi']
    ZODI_LEVELS = {  #Normalized to Johnson V
        'None': None,
        'Low': 23.3,
        'Average': 22.7,
        'High': 22.1
}

    #TODO: handle the posfile in the same way as the other data files
    #  except that some variable expansion has to be done, unlike for
    #  the other files for which that is handled by pysynphot
    #posfile = MODELS['zodipos']

    def __init__(self, formdata=None):
        #Then call the parent
        SkyElement.__init__(self, formdata)
        #This tells us which type of calculation we're going to do
        self.ztype = self.fd['ZodiSpec']
        self._set_norm()

    def make_synexpr(self):
        if self.norm is None:
            self.expr = None
        else:
            if self.ztype in ('ZodiStandard', 'ZodiMag', 'ZodiCoords'):
                self.expr = "rn(%s,band(johnson,v),%f,vegamag)"%(self.specfile,self.norm)
            elif self.ztype in ('ZodiMult'):
                self.expr = "%f*%s"%(self.norm,self.specfile)

    def _set_norm(self):
        #There are many ways to specify how the Zodiacal Light spectrum should be normalized
        if self.ztype == 'ZodiStandard':
            level = self.fd[self.ztype]
            self.norm = self.ZODI_LEVELS[level]
            self.descr = level

        elif self.ztype == 'ZodiMag':
            self.norm = self.fd[self.ztype]
            self.descr = 'User Defined (normalized to %f V)'%self.norm

        elif self.ztype == 'ZodiMult':
            self.norm = self.fd[self.ztype]
            self.descr = 'User Defined (scaled by %f)'%self.norm

        elif self.ztype == 'ZodiCoords':
            zra, zdec = self.fd['ZodiRA'], self.fd['ZodiDec']
            ra, dec = angles.radec_to_decimal( zra, zdec )

            if self.fd['ZodiSunAttribute'] == 'Date':
                year, month, day = self.fd['ZodiDate'], self.fd['ZodiMonth'], self.fd['ZodiDay']
                date = datetime.date(year,month,day)
                helong, helat, sunangle = angles.radec_date_to_helonglat(ra, dec, date)
                self.descr = ('%s %s %s/%s/%s (%f deg.)' % (zra,zdec,month,day,year,sunangle))

            elif self.fd['ZodiSunAttribute'] == 'SunAngle':
                zsun = self.fd['ZodiSun']
                # Table is defined from 0-180, but users can put in anything
                sunangle = angles.limit_longitude(float(zsun))
                self.descr = ('%s %s Sun Angle = %f deg.' % (zra, zdec, zsun))
                helong, helat = angles.radec_sunangle_to_helonglat(ra, dec, sunangle)

            else:
                raise ValueError('Unknown sun specification %s'%self.fd['ZodiSunAttribute'])
            try:
                self.norm = angles.zodi_lookup(helong, helat)
            except RangeError as e:
                e.addinfo("These values are derived from the RA, Dec, and Date or Helio-ecliptic Longitude. You must either select different values for these quantities, or a different normalization method for the zodiacal light.")
                raise

        else:
            raise ValueError('Unknown normalization type %s' % self.ztype)


class Sky(object):
    """Class to compute the sky background including all relevant SkyElements.
    Presently supports these elements:
        Airglow, Earthshine, HeliumEmiss, and ZodiacalLight.
    Not all the elements are required. Elements will be picked out based on
    the contents of the dictionary that is passed in."""
    def __init__(self, formdata=None):
        self.airglow=None
        self.earthshine=None
        self.helium=None
        self.zodi=None
        self.descr=dict()
        self.expr='Not yet set'

        if formdata:
            try:
                self.airglow=Airglow(formdata)
                self.airglow.make_synexpr()
                self.descr['airglow_str']=self.airglow.descr
            except KeyError as e:
                pass

            try:
                self.earthshine=Earthshine(formdata)
                self.earthshine.make_synexpr()
                self.descr['earthshine_str']=self.earthshine.descr
            except KeyError as e:
                pass
            try:
                self.helium=HeliumEmiss(formdata)
                self.helium.make_synexpr()
                self.descr['helium_str']=self.helium.descr
            except KeyError as e:
                pass

            try:
                self.zodi= ZodiacalLight(formdata)
                self.zodi.make_synexpr()
                self.descr['zodi_str']=self.zodi.descr
            except KeyError as e:
                pass

    def __str__(self):
        return str(self.descr)

    def make_synexpr(self):
        elements = []
        for x in (self.airglow, self.earthshine, self.helium, self.zodi):
            if x and x.expr:
                elements.append(x.expr)
        if len(elements) > 0:
            self.expr = ' + '.join(elements)
        else:
            #Provide a valid expression with a value of zero
            self.expr = 'unit(0.0,flam)'
