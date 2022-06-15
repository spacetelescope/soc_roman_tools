'''Map the sky background form parameters into a form that the etc_engine can use:
(ie, a synphot expression)
keep this file in sync with 
https://github.com/spacetelescope/pyetc/blob/master/pyetc/etc_web/etc/sky_angles.py
'''

from __future__ import division
import os
import datetime
import numpy as np

from pandeia.engine.helpers.background.hst import util
from pandeia.engine.custom_exceptions import BackgroundError, RangeError

DAYS_IN_YEAR = 365.25
OBLIQUITY = 23.43 / 180.0 * np.pi


def radec_to_elonglat(ra, dec):
    '''convert ra, dec to ecliptic long,long (all in decimal deg)'''

    # degrees to radians
    ra = ra / 180.0 * np.pi
    dec = dec / 180.0 * np.pi

    # ecliptic latitude
    sin_lat = np.sin(dec) * np.cos(OBLIQUITY) - np.cos(dec) * np.sin(OBLIQUITY) * np.sin(ra)
    lat = np.arcsin(sin_lat)

    # ecliptic longitude
    tan_long = (np.tan(dec) * np.sin(OBLIQUITY) + np.cos(OBLIQUITY) * np.sin(ra)) / np.cos(ra)
    long = np.arctan(tan_long)

    # convert longitude to appropriate quadrant
    cos_long = np.cos(dec) * np.cos(ra) / np.cos(lat)
    sin_long= (np.sin(dec) * np.sin(OBLIQUITY) + (np.cos(dec) * np.sin(ra) * np.cos(OBLIQUITY))) / np.cos(lat)
    long_cos = np.arccos(cos_long)
    long_sin = np.arcsin(sin_long)
    if long < 0 and long_sin > 0 and long_cos > 0:
        long += np.pi
    elif long > 0 and long_sin < 0 and long_cos > 0:
        long += np.pi
    elif long < 0 and long_sin < 0 and long_cos > 0:
        long += np.pi * 2.

    # radians to degrees,
    lat  = lat * 180.0 / np.pi
    long = long * 180.0 / np.pi

    return long, lat


def radec_date_to_helonglat(ra, dec, date):
    '''compute the heliocentric ecliptic longitude/latitude given ra, dec and obs date

    This is not a high precision routine!!  (since zodical variance with angle is slow)
    Input and output angles are in degrees, date is a datetime object. Angle
    computations done in radians, however.
    '''
    # determine the longitude of the vernal equinox given the date
    # calculate fraction of year separating the dates
    veq = datetime.date(date.year,3,20)
    delta = (date - veq).days
    if delta > DAYS_IN_YEAR:
        delta -= DAYS_IN_YEAR
    sunangle = 360.*delta/DAYS_IN_YEAR
    elong, elat = radec_to_elonglat(ra, dec)

    # Don't limit the ranges of elong and sunangle (using
    # limit_longitude) until ready to look up in the table.
    return elong-sunangle, elat, sunangle

def limit_longitude(theta):
    theta = np.abs(theta)
    if theta > 180.0:
        return 360.0 - theta
    else:
        return theta

def limit_latitude(theta):
    theta = np.abs(theta)
    if theta > 90.0:
        return 180.0 - theta
    else:
        return theta

def radec_sunangle_to_helonglat(ra, dec, sunangle):
    '''compute the heliocentric ecliptic longitude/latitude given ra, dec and longitudnal
    angle to sun

    This is not a high precision routine!!  (since zodical variance with angle is slow)
    Input and output angles are in degrees, date is a datetime object. Angle
    computations done in radians, however.
    '''
    elong, elat = radec_to_elonglat(ra, dec)
    # Don't limit the ranges of elong and sunangle (using
    # limit_longitude) until ready to look up in the table.
    return sunangle, elat

def radec_to_decimal(in_ra, in_dec):
    """Convert in_ra, in_dec from input values to floating-point degrees; returns ra, dec

    For each, if we can convert it to a float, use the result.  Otherwise,
    parse HH:MM:SS or DD:MM:SS to convert it
    """

    try:
        ra = float(in_ra)
    except ValueError:
        # Then it must be hours:minutes:seconds
        l = in_ra.split(':')

        # convert to hours
        ra = int(l[0])
        if len(l) > 1 :
            ra += float(l[1])/60.
        if len(l) > 2 :
            ra += float(l[2])/3600.
        #
        # convert hours to degrees
        ra *= 15.

    try:
        dec = float(in_dec)
    except ValueError:
        # then it must be degrees:minutes:seconds
        l = in_dec.split(':')

        # get unsigned degrees
        dec = np.abs(int(l[0]))
        if len(l) > 1 :
            dec += float(l[1])/60.
        if len(l) > 2 :
            dec += float(l[2])/3600.

        # apply the sign to it.
        if in_dec.strip().startswith('-') :
            dec = -dec

    return ra, dec

def bilinear_interp(xref, yref, arr, x, y):
    '''perform bilinear interpolation on 2-d array.

    xref: contains x values corresponding to first dim index in array.
    yref: contains y values corresponding to second dim index in array.
    arr: array of values to be interpolated.
    x, y: x and y values for interpolation, scalars or arrays permitted.

    Poor error handling at the moment...
    '''

    xind = np.searchsorted(xref, x)
    yind = np.searchsorted(yref, y)

    # Bound the indices to 0 .. N-1
    xind0 = np.maximum(0, xind - 1)
    yind0 = np.maximum(0, yind - 1)
    xind1 = np.minimum(len(xref) - 1, xind)
    yind1 = np.minimum(len(yref) - 1, yind)

    # The denominator may be zero if outside of the bounds of the
    # image
    xden = (xref[xind1] - xref[xind0])
    xden = np.where(xden == 0.0, 1.0, xden)
    yden = (yref[yind1] - yref[yind0])
    yden = np.where(yden == 0.0, 1.0, yden)

    xfract = (x - xref[xind0]) / xden
    xfracti = 1.0 - xfract
    yfract = (y - yref[yind0]) / yden
    yfracti = 1.0 - yfract

    val = (arr[xind0,yind0] * xfracti * yfracti +
           arr[xind1,yind0] * xfract * yfracti +
           arr[xind0,yind1] * xfracti * yfract +
           arr[xind1,yind1] * xfract * yfract)
    return val

def zodi_lookup(helong, helat, posfile=None):
    '''look up zodical magnitude given ecliptic latitude and longitude in deg.'''
    #This lookup table presently lives with the source code
    if posfile is None:
        posfile = os.path.join(os.path.dirname(__file__),
                               'zodi.dat')
    edict = util.read_dict(posfile)

    # The table is symmetric for 0 <= long <= 180 and 0 <= lat <= 90,
    # so limit the lookup values to those ranges.
    helong = limit_longitude(helong)
    helat = limit_latitude(helat)

    table = np.array(edict['table'])
    marr = table[1:,1:] #table contents
    xref = table[1:,0]  #header (index) row
    yref = table[0,1:]  #header (index) column
    # look up value
    zmag = bilinear_interp(xref, yref, marr, np.abs(helong), np.abs(helat))

    if np.isnan(zmag):
        raise RangeError('Helio-ecliptic long and ecliptic lat (%f, %f) outside valid range of zodiacal light table.'
                         %(helong,helat))
    return zmag
