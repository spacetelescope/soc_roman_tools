# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import
from astropy import units as u
import synphot as syn


# These are the wavelength and flux units that are used internally for all calculations.
PANDEIA_WAVEUNITS = u.micron
PANDEIA_FLUXUNITS = u.mJy

UNIT_MAP = {"flam": syn.units.FLAM, "fnu": syn.units.FNU, "abmag": u.ABmag, 
            "vegamag": syn.units.VEGAMAG, 'njy': u.nJy, "ujy": u.uJy, "mjy": u.mJy,
            "jy": u.Jy, "MJy": u.MJy, "microns": u.micron, "micron": u.micron, 
            "angstrom": u.AA, "nm": u.nm, "um": u.um, "mm": u.mm, "m": u.m}

# Absolute and relative tolerance values to use to check if values are close enough to be considered equal.
# See https://docs.scipy.org/doc/numpy/reference/generated/numpy.isclose.html for details on how they're applied
# in comparisons.
PANDEIA_ATOL = 2.0e-8
PANDEIA_RTOL = 1.0e-6

# Pandeia subsamples lines in order to accurately reproduce thin spectral lines in a spectrum. This subsampling must
# extend over a large enough range to contain most of the wings and accurately reproduce the line.
#
# spectral_line_window determines the size of the wavelength window over which lines are subsampled. Values above 15
#  (decided empirically, see Issue #3292) are necessary to contain the wings of the line.
# spectral_line_subsample determines how much the spectrum is subsampled within the window to reproduce the shape of
#  the profile.
SPECTRAL_LINE_WINDOW = 20
SPECTRAL_LINE_SUBSAMPLE = 3

# In practice, a value of 200 samples within an imaging configuration's wavelength range (i.e. filter bandpass) should
# be more than enough to produce an accurate image. Note that because we use synphot to resample, even the flux of
# narrow lines is conserved.
SPECTRAL_MAX_SAMPLES = 200

# Pandeia's spectrum processing can use a considerable amount of memory for large fields (like SOSS). To handle this,
# we divide the spectra into chunks of wavelength planes (if necessary) at the cost of a few seconds of run time. There
# may be a more elegant way to handle this...
SPECTRAL_MAX_CHUNK = 1000

# Pandeia's FOV size (for the spatial axes of the scene cubes) is the largest of either the scene size, the PSF size,
# or the configured instrument's default FOV size. We need to add a small buffer to the scene size to adequately contain
# the full scene after sources have been represented by PSFs
SCENE_FOV_PIXBUFFER = 20

# The ConvolvedSceneCube function uses astropy.convolve_fft to produce a convolved scene cube. This cube cannot be
# created if any of the pixels are too bright - even getting close to the maximum floating point number that can be
# stored in a numpy array. Empirical studies show that values greater than 1e33 still cause astropy.convolve_fft to fail
# We have set this to 1e32 to be safe.
SCENE_MAX_PIXELVAL = 1e32

# The profile generator needs to find pixels whose locations are 0, but due to floating point math, the values won't
# necessarily be exactly zero. Set a tiny value accordingly.
POS_EPSILON = 1e-8

# The profile generator needs to oversample the central pixel of high-sersic-index sersic profiles. They can be
# incredibly peaky, such that the amount of flux would be drastically overstated if we were to simpy set the central
# pixel to the value at the exact center. Instead, we need to oversample that pixel and compute a more accurate flux for
# the central pixel of the sersic. The process by which the number 101 was arrived at can be found here:
# Code PR #3892#issuecomment-397713379
SERSIC_OVERSAMPLE = 101

# The default dispersion is meant for use when attempting to create backgrounds for order 0 of slitless and multiorder
# modes, because order 0 itself is not dispersed but our algorithms demand it.
DEFAULT_DISPERSION = 0.01

# Minimum to clip values to, for avoidance of divsion-by-zero errors in exposure timing, profile generation, etc.
MIN_CLIP = 1e-10
