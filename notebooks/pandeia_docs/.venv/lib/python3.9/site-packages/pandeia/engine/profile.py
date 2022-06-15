# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import numpy as np
import scipy.special as sp
import astropy.units as u
from scipy.integrate import dblquad

from .config import DefaultConfig
from .coords import Grid
from .constants import POS_EPSILON, SERSIC_OVERSAMPLE, MIN_CLIP
from .utils import recursive_subclasses
from .custom_exceptions import EngineInputError


def ProfileFactory(config={}, **kwargs):

    """
    Function to take profile configuration data and build/return a configured instance of
    the appropriate distribution generator.

    Parameters
    ----------
    webapp: bool
        Switch to toggle strict API checking
    config: dict
        Configuration data in engine API dict format
    **kwargs: keyword/value pairs
        Additional configuration data
    """
    #all_config = merge_data(config, dict(**kwargs))
    all_config = config
    if 'geometry' in all_config.shape:
        method = all_config.shape['geometry'].replace('_', '')
    else:
        msg = "Must specify type of profile to create."
        raise EngineInputError(value=msg)

    # we need to run the API checks here after merging in the normalization type-specific defaults
    # if webapp:
    #    self._api_checks(all_config)

    types = recursive_subclasses(Distribution)
    methods = [t.__name__.lower().replace('distribution', '') for t in types]
    type_map = dict(list(zip(methods, types)))

    if method not in methods:
        msg = "Unsupported or not yet implemented profile type: %s" % method
        raise EngineInputError(value=msg)
    else:
        cls = type_map[method](config,**kwargs)
        return cls


class Distribution(DefaultConfig):
    """
    Abstract class for source profiles.
    Normally, when defining undefined abstract properties, we use NotImplementedError, with an eye toward implementing
    that feature later. However, in this case, when the normalization methods are not implemented it is because there
    is no robust sensible definition for them, and we will NEVER implement them. Therefore, the abstract classes
    produce errors to that effect.

    Parameters
    ----------
    src: dict
        Details of the source geometry to generate, from the Source block of the input files.

    Attributes
    ----------
    xoff: float
        Offset in X direction
    yoff: float
        Offset in Y direction
    pa: float
        Position angle of major axis
    xrot: np.ndarray
        rotated and shifted indicies, X axis values
    yrot: np.ndarray
        rotated and shifted indicies, Y axis values
    pix_area_sqarcsec: float
        area of a detector pixel in square arcseconds
    prof: np.ndarray
        raw 2D profile in detector units
    norm_prof: np.ndarray
        normalized 2D profile in detector units

    Methods
    -------
    normalized()
        return the normalized 2D profile
    raw()
        return the raw 2D profile

    """
    def __init__(self,src):
        self.pa=src.position['orientation']
        self.xoff=src.position['x_offset']
        self.yoff=src.position['y_offset']
        self.grid = src.grid
        yrot, xrot = self.grid.shift_rotate(self.yoff, self.xoff, self.pa)
        self.xrot = xrot
        self.yrot = yrot
        self.pix_area_sqarcsec = src.grid.xsamp * src.grid.ysamp
        self.norm_method = src.shape['norm_method']
        self.geometry = src.shape['geometry']

    def pixelscale(self):
        if self.surf_area_units in ['sr']:
            arcsec2 = u.arcsec * u.arcsec
            normfactor = self.pix_area_sqarcsec / u.sr.to(arcsec2)  # convert area in steradians to area in pixels
        elif self.surf_area_units in ['arcsec^2', None]: # 'None' should be an option because integrated flux
                                                        # shouldn't have units (internally, the grid is arcsec)
            normfactor = self.pix_area_sqarcsec
        else:
            msg = "Unsupported surface area unit: %s" % self.surf_area_units
            raise EngineInputError(value=msg)

        return normfactor
    def normalize(self):
        options = {'integ_infinity': self.integrate_infinity, 'surf_scale': self.surface_scale,
                       'surf_center': self.surface_center}
        options[self.norm_method]()

    def integrate_infinity(self):
        msg = "Normalization method {0:} not supported for profile type {1:}".format(self.norm_method, self.geometry)
        raise EngineInputError(value=msg)

    def surface_scale(self):
        msg = "Normalization method {0:} not supported for profile type {1:}".format(self.norm_method, self.geometry)
        raise EngineInputError(value=msg)

    def surface_center(self):
        msg = "Normalization method {0:} not supported for profile type {1:}".format(self.norm_method, self.geometry)
        raise EngineInputError(value=msg)

    def normalized(self):
        return self.norm_prof

    def raw(self):
        return self.prof


class PointDistribution(Distribution):
    """
    Use pandeia.engine.coords.Grid.point_source() to generate a point source with
    subpixel positioning.

    Parameters
    ----------
    src: dict
        Details of the source geometry to generate, from the Source block of the input files.

    Attributes
    ----------
    xoff: float
        Offset in X direction
    yoff: float
        Offset in Y direction
    prof: np.ndarray
        raw 2D profile in detector units
    norm_prof: np.ndarray
        normalized 2D profile in detector units

    """
    def __init__(self, src):
        self.xoff = src.position['x_offset']
        self.yoff = src.position['y_offset']

        # This function both creates and fills the point source profile. It should
        # eventually be refactored to work the same way as other sources
        self.prof = src.grid.point_source(xoff=self.xoff, yoff=self.yoff)

        self.norm_prof = self.prof

class SersicDistribution(Distribution):
    """
    Create a 2-dimensional elliptical source on the current grid. The intensity
    profile is described by a Sersic profile, I(r) = I(0) * exp(-b*((r/r_scale)**(1/n))-1),
    where r_scale is the effective radius such that I(r) = 0.5*I(0) and n is the Sersic index.
    The ellipticity is governed by specifying major and minor axis scale lengths
    separately.


    Parameters
    ----------
    src: dict
        Details of the source geometry to generate, from the Source block of the input files.

    Attributes
    ----------
    minor: float
        Minor axis scale length
    major: float
        Major axis scale length
    pa: float
        Position angle in degrees of major axis measured positive in +X direction
    xoff: float
        Offset in X direction
    yoff: float
        Offset in Y direction
    sersic_index: float
        Sersic profile shape parameter. 0.5 => gaussian, 1.0 => exponential, 4.0 => de Vaucouleurs
    prof: np.ndarray
        raw 2D profile in detector units
    norm_prof: np.ndarray
        normalized 2D profile in detector units

    Methods
    -------
    _integrate_infinity():
        We use an analytical function to represent the integration of the sersic function
        from -Inf to +Inf in both axes. This will account for flux that falls
        outside of the FOV.
    _surface_center():
        Normalization is to the surface brightness of the center of the profile. This
        formulation is already normalized at that point.
    _surface_scale():
        Normalization is to the surface brightness of a point at the e-folding radius of
        this profile, where intensity is 1/e of the central intensity.
    sersic_generator():
        Handle subpixel-positioning and proper modeling of the central pixel (which can be too
        bright when sersic_index > 3)
    """
    def __init__(self,src):
        Distribution.__init__(self,src)
        self.norm_method = src.shape['norm_method']
        self.major=src.shape['major']
        self.minor=src.shape['minor']
        self.sersic_index=src.shape['sersic_index']
        self.surf_area_units = src.shape['surf_area_units']
        # the actual value of b. Formula taken from astropy's sersic2d shape.
        self.b = sp.gammaincinv(2*self.sersic_index,0.5)

        self.sersic_generator()

        self.normalize()

    def sersic_generator(self):
        '''
        Generic function to create a properly normalized sersic, taking into account:
         * Subpixel positioning (off-screen and on) by replicating some of the grid.point_source functionality
         * Realistic peak counts, by creating a highly oversampled grid covering just the central pixel

        NB: Subclasses should override sersic_func to specify the appropriate equation for the desired distribution.
        This method will rely on the subclass's sersic_func method.

        '''
        # oversample the central pixel to compute an actual flux for what may be an extremely tiny peak
        center_grid = Grid(self.grid.xsamp/(SERSIC_OVERSAMPLE-1), self.grid.ysamp/(SERSIC_OVERSAMPLE-1),
                           SERSIC_OVERSAMPLE, SERSIC_OVERSAMPLE)
        yrot,xrot = center_grid.shift_rotate(0, 0, self.pa)
        center_pixel = self.sersic_func(yrot,xrot,self.major,self.minor,self.sersic_index)
        # this should only correct for the oversampling compared to the usual grid
        self.center_flux = np.sum(center_pixel) / (SERSIC_OVERSAMPLE**2)

        # self.prof is a container to hold the actual profile
        self.prof = np.zeros_like(self.grid.x)
        # find the pixel-snapped locations of the lower left of the four pixels, and how far off the actual center it is
        # note that this doesn't actually have to be on the detector grid, this is just to properly split up the flux.
        lx = (self.xoff // self.grid.xsamp) * self.grid.xsamp
        rx = (self.xoff % self.grid.xsamp) * self.grid.xsamp
        ly = (self.yoff // self.grid.ysamp) * self.grid.ysamp
        ry = (self.yoff % self.grid.ysamp) * self.grid.ysamp

        # make a list of pixel positions and the amount of flux distributed to each.
        pixels = [[lx, ly, (1 - rx) * (1 - ry)], [lx, ly + self.grid.ysamp, (1 - rx) * ry],
                      [lx + self.grid.xsamp, ly, rx * (1 - ry)],
                      [lx + self.grid.xsamp, ly + self.grid.ysamp, rx * ry]]
        # create profiles centered in each of the four pixels, and scale them to the amount of flux
        for tx,ty,frac in pixels:
            # skip pixel if it adds no flux. In theory, the flux will be split between 1, 2, or 4 pixels.
            # if it's exactly centered on a pixel (say, at 0,0) the flux will be entirely on one pixel.
            # if it's offset slightly in X or Y, it will be split across two pixels
            # if it's offset slightly in both X and Y, it'll be split among four pixels.
            if frac != 0:
                # xrot and yrot contain the coordinates of each detector pixel center, in the coordinate system of the
                # offset and centered profile.
                yrot, xrot = self.grid.shift_rotate(ty, tx, self.pa)
                partial = self.sersic_func(yrot, xrot, self.major, self.minor, self.sersic_index)
                # the zero-value may not be zero, but something exceedingly small
                yzero = np.min(np.abs(yrot))
                xzero = np.min(np.abs(xrot))
                # if the central pixel is actually on the detector grid (and the sersic_index is greater than 0.5), we
                # need to find it and replace its flux with the corrected flux for the central pixel
                if np.abs(xzero) < POS_EPSILON * self.grid.xsamp and np.abs(yzero) < POS_EPSILON * self.grid.ysamp:
                    xvals = np.where(np.abs(xrot.flatten()) < POS_EPSILON * self.grid.xsamp)[0]
                    yvals = np.where(np.abs(yrot.flatten()) < POS_EPSILON * self.grid.ysamp)[0]
                    xy = np.intersect1d(xvals,yvals)
                    xval = xy // self.grid.shape[0]
                    yval = xy % self.grid.shape[0]
                    partial[xval, yval] = self.center_flux
                self.prof += partial * frac

    def integrate_infinity(self):
        # integrate the Sersic profile to get the total flux for normalization, including flux outside the FOV
        # http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
        integral = self.major * self.minor * 2 * np.pi * self.sersic_index * np.exp(self.b)/(self.b**(2*self.sersic_index))* sp.gamma(2 * self.sersic_index)
        self.norm_prof = self.prof / integral * self.pixelscale()

    def surface_scale(self):
        # for the base sersic profile, the flux at R=Re is 1. e^(-b* (1-1))
        self.norm_prof = self.prof * self.pixelscale()

    def surface_center(self):
        # for the base sersic profile, the flux at R=0 is e**b, so we must multiply by e**-b
        self.norm_prof = self.prof * self.pixelscale() * np.e**(-1*self.b)

    def sersic_func(self, y, x, major, minor, index):
        """
        Implement the Sersic intensity profile as a function to actually fill the grid with the profile.
        This is Equation 1 of Graham & Driver (2005) 2005PASA...22..118G

        This distribution function is also replicated in client/js/scenepage.js and should be kept in sync with any
        changes made here.

        Parameters
        ----------
        y: float or numpy.ndarray
            Y values for evaluating function
        x: float or numpy.ndarray
            X values for evaluating function
        major: float
            Major axis scale length
        minor: float
            Minor axis scale length
        index: float
            Sersic index

        Returns
        -------
        profile: float or numpy.ndarray
            Float or array containing evaluated Sersic profile
        """
        dist = np.sqrt((x / major)**2.0 + (y / minor)**2.0)
        # This is Equation 1 of Graham & Driver (2005) 2005PASA...22..118G
        profile = np.exp( -self.b * (dist**(1.0 / index) - 1) )

        return profile


class SersicScaleDistribution(SersicDistribution):
    """
    Create a 2-dimensional elliptical source on the current grid. The intensity
    profile is described by a Sersic profile, I(r) = I(0) * exp(-(r/r_scale)**(1/n)),
    where r_scale is the scale length where I(r) = I(0)/e and n is the Sersic index.
    The ellipticity is governed by specifying major and minor axis scale lengths
    separately.
    This is Equation 14 of Graham & Driver (2005) 2005PASA...22..118G


    This was the original implementation of the Sersic profile in pandeia, up through version 1.2.2. However, based on
    community feedback, it was replaced by the above SersicDistribution class. This class was retained but renamed to
    permit backwards compatibility in the web app.


    Parameters
    ----------
    src: dict
        Details of the source geometry to generate, from the Source block of the input files.

    Attributes
    ----------
    minor: float
        Minor axis scale length
    major: float
        Major axis scale length
    pa: float
        Position angle in degrees of major axis measured positive in +X direction
    xoff: float
        Offset in X direction
    yoff: float
        Offset in Y direction
    sersic_index: float
        Sersic profile shape parameter. 0.5 => gaussian, 1.0 => exponential, 4.0 => de Vaucouleurs
    prof: np.ndarray
        raw 2D profile in detector units
    norm_prof: np.ndarray
        normalized 2D profile in detector units

    Methods
    -------
    _integrate_infinity():
        We use an analytical function to represent the integration of the sersic function
        from -Inf to +Inf in both axes. This will account for flux that falls
        outside of the FOV.
    _surface_center():
        Normalization is to the surface brightness of the center of the profile. This
        formulation is already normalized at that point.
    _surface_scale():
        Normalization is to the surface brightness of a point at the e-folding radius of
        this profile, where intensity is 1/e of the central intensity.
    """

    def __init__(self, src):
        Distribution.__init__(self, src)
        self.norm_method = src.shape['norm_method']
        self.major = src.shape['major']
        self.minor = src.shape['minor']
        self.sersic_index = src.shape['sersic_index']
        self.surf_area_units = src.shape['surf_area_units']

        # The sersic_generator is inherited from the parent SersicDistribution class, but will correctly use this
        # class's own sersic_func(), which contains the appropriate equation.
        self.sersic_generator()

        self.normalize()

    def integrate_infinity(self):
        # integrate the Sersic profile to get the total flux for normalization, including flux outside the FOV
        # http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
        integral = self.major * self.minor * 2 * np.pi * self.sersic_index * sp.gamma(2*self.sersic_index)
        self.norm_prof = self.prof / integral * self.pixelscale()

    def surface_scale(self):
        self.norm_prof = self.prof * np.e * self.pixelscale()

    def surface_center(self):
        self.norm_prof = self.prof * self.pixelscale()

    def sersic_func(self, y, x, major, minor, index):
        """
        Implement the Sersic_scale intensity profile as a function to actually fill the grid with the profile.

        This distribution function is also replicated in client/js/scenepage.js and should be kept in sync with any
        changes made here.

        Parameters
        ----------
        y: float or numpy.ndarray
            Y values for evaluating function
        x: float or numpy.ndarray
            X values for evaluating function
        major: float
            Major axis scale length
        minor: float
            Minor axis scale length
        index: float
            Sersic index (traditionally, a number between 0.5 (gaussian) and 10)

        Returns
        -------
        profile: float or numpy.ndarray
            Float or array containing evaluated Sersic profile
        """
        dist = np.sqrt((x / major) ** 2.0 + (y / minor) ** 2.0)
        # This is Equation 14 of Graham & Driver (2005) 2005PASA...22..118G
        profile = np.exp(-dist ** (1.0 / index))
        return profile


class Gaussian2dDistribution(SersicScaleDistribution):
    """
    The Gaussian 2D profile is a special case of the Sersic_scale profile, where sersic index = 0.5.

    Parameters
    ----------
    src: dict
        Details of the source geometry to generate, from the Source block of the input files.

    Attributes
    ----------
    minor: float
        Minor axis scale length
    major: float
        Major axis scale length
    pa: float
        Position angle in degrees of major axis measured positive in +X direction
    xoff: float
        Offset in X direction
    yoff: float
        Offset in Y direction
    prof: np.ndarray
        raw 2D profile in detector units
    norm_prof: np.ndarray
        normalized 2D profile in detector units

    """
    def __init__(self,src):
        Distribution.__init__(self,src)
        self.norm_method = src.shape['norm_method']
        self.major=src.shape['major']*np.sqrt(2.0) # to match the usual definition of a gaussian
        self.minor=src.shape['minor']*np.sqrt(2.0)
        self.sersic_index=0.5
        self.surf_area_units = src.shape['surf_area_units']

        # The sersic_generator is inherited from the parent SersicDistribution class, but will correctly use
        # this class's (SersicScaleDistribution) own sersic_func() which contains the appropriate equation for
        # a gaussian2d function.
        self.sersic_generator()

        self.normalize()


class FlatDistribution(Distribution):
    """
    Implement a source with a constant surface brightness as an ellipse, i.e. a tilted disc.

    Parameters
    ----------
    src: dict
        Details of the source geometry to generate, from the Source block of the input files.

    Attributes
    ----------
    minor: float
        Minor axis scale length
    major: float
        Major axis scale length
    pa: float
        Position angle in degrees of major axis measured positive in +X direction
    xoff: float
        Offset in X direction
    yoff: float
        Offset in Y direction
    prof: np.ndarray
        raw 2D profile in detector units
    norm_prof: np.ndarray
        normalized 2D profile in detector units

    Methods
    -------
    integrate_infinity():
        Normalization to flux at infinity means total flux, and that if all pixels in the
        profile are summed up, they should sum to 1.
    surface_center():
        Normalization to the center is, in this case, a no-op, because the initial
        definition of the disk sets the center (and everywhere else) to 1

    """
    def __init__(self,src):
        Distribution.__init__(self,src)
        self.norm_method = src.shape['norm_method']
        self.major=src.shape['major']
        self.minor=src.shape['minor']
        self.surf_area_units = src.shape['surf_area_units']

        self.prof = self.grid.elliptical_mask(self.major, self.minor, pa=self.pa, xoff=self.xoff, yoff=self.yoff)

        self.normalize()

    def integrate_infinity(self):
        self.norm_prof = self.prof / (np.pi * self.major * self.minor) * self.pixelscale()

    def surface_center(self):
        self.norm_prof = self.prof * self.pixelscale()


class PowerDistribution(Distribution):
    """
    Create a 2-dimensional power law distribution on the current grid. I(r) = I(0) * r**-k
    There is only one configurable parameter, k, which must be positive.
    It has the same basic shape as a power law function, apart from a flat core within a radius r_core, which sidesteps
    the singularity at r=0.

    This function cannot be integrated to infinity unless k > 2. Because there is a scientific desire to have powers
    between 0 and 2, we have opted to not allow integration to infinity at all.

    Parameters
    ----------
    src: dict
        Details of the source geometry to generate, from the Source block of the input files.

    Attributes
    ----------
    power_index: float
        The slope of the power law index. Must be positive.
    r_core: float
        The radius of the flat core of the power law profile (arcsec)
    prof: np.ndarray
        raw 2D profile in detector units
    norm_prof: np.ndarray
        normalized 2D profile in detector units

    Methods
    -------
    surface_center():
        Normalization is to the surface brightness of the center of the profile. This
        formulation is already normalized at that point.
    """

    def __init__(self, src):
        Distribution.__init__(self, src)
        self.power_index = src.shape['power_index']
        self.surf_area_units = src.shape['surf_area_units']
        self.r_core = src.shape['r_core']
        if self.power_index <= 0:
            raise ValueError('Power Law Index must be positive, not {}'.format(self.power_index))

        self.prof = self.exponential_func(self.yrot, self.xrot, self.r_core, self.power_index)

        self.normalize()

    #def integrate_infinity(self):
    #    # integrate the Power Law profile to get the total flux for normalization, including flux outside the FOV
    #    # This is two integrations: One over the central circle; one for the exponential profile outside of that.
    #    integral = np.pi * self.r_core**2 + 2* np.pi * self.r_core**2/(self.power_index - 2)
    #    self.norm_prof = self.prof / integral * self.pixelscale()

    def surface_center(self):
        # the profile comes pre-normalized to 1 at the center.
        self.norm_prof = self.prof * self.pixelscale()

    def exponential_func(self, y, x, r_core, index):
        """
        Implement an exponential function. This isn't a pure exponential, because the real one is infinite at dist=0
        but this shares the same basic properties.

        Any changes to this profile will need to be made to ui/client/js/scenepage.js as well.

        Parameters
        ----------
        y: float or numpy.ndarray
            Y values for evaluating function
        x: float or numpy.ndarray
            X values for evaluating function
        r_core: float
            Radius of the flat core
        index: float
            power law index
        Returns
        -------
        profile: float or numpy.ndarray
            Float or array containing evaluated Sersic profile
        """
        dist = np.sqrt((x/r_core)**2.0 + (y/r_core)**2.0)
        profile = (dist.clip(MIN_CLIP, np.max(dist)))**(-1*index)
        # flatten the central portion. Everything within the core radius is set to 1.
        profile[np.where(dist <= 1.0)] = 1.0

        return profile
