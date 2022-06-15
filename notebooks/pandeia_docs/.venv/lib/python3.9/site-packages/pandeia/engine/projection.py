# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import numpy as np

from .coords import Grid, IrregularGrid
from .utils import recursive_subclasses
from .custom_exceptions import EngineOutputError, EngineInputError

def ProjectionFactory(projection_type):
    """
    Function to take projection type and return the appropriate projection-dependent functions.

    Parameters
    ----------
    projection_type: string
        String with four valid values: "Image", "Spec", "Slitless", "Multiorder".
    """

    types = recursive_subclasses(Projection)
    projections = [t.__name__.lower().replace('projection', '') for t in types]
    type_map = dict(list(zip(projections, types)))

    if projection_type not in projections:
        msg = "Unsupported or not yet implemented projection: %s" % projection_type
        raise EngineInputError(value=msg)
    else:
        cls = type_map[projection_type]()
        return cls


class Projection(object):
    """
    Abstract class with projection-type dependent code. Currently only includes 
    projection-dependent aspects of the Report class, but should eventually be expanded to 
    include projection-dependent parts of signal.py, noise.py, and possibly instrument.py

    Note that none of these correspond to the actual 2D results of the IFU projection type, 
    which technically isn't a projection as it reports the full scene cubes.

    It also does not include Coronagraphy, which (uniquely) has no background to subtract 
    from its _2d results.

    """

    def __init__(self):
        pass

    def _2d(self, extracted, saturation, signal):

        self.signal = signal

        # Signal and noise in 2D. Get from extracted products
        s = extracted['detector_signal']
        n = extracted['detector_noise']
        t = extracted['detector_saturation']

        # This is the background rate in each pixel without sources
        self.bg_pix = signal.rate_plus_bg - signal.rate

        # The noise is undefined if the pixel is fully saturated.
        # The engine currently returns noise=0 if the pixel has full saturation, which is
        # confusing since the noise is not 0, but rather undetermined or infinite. Setting
        # the noise to NaN ensures that the S/N of saturated pixels are NaNs.
        n[t == 2] = np.nan
        self.detector_sn = (s - self.bg_pix) / n
        self.detector_signal = s + n * np.random.randn(n.shape[0], n.shape[1])
        self.saturation = t
        self.ngroups_map = extracted['detector_ngroups']


class ImageProjection(Projection):

    def _2d(self, extracted, saturation, signal):
        Projection._2d(self,extracted,saturation,signal)

        return {'detector': self.detector_signal, 'snr': self.detector_sn, 'saturation': self.saturation, \
                'ngroups_map': self.ngroups_map}, self.bg_pix

    def _project(self, signal, extracted):

        return signal.grid, signal.wave_pix, signal.pixgrid_list[0]


class SpecProjection(Projection):

    def _2d(self, extracted, saturation, signal):
        Projection._2d(self,extracted,saturation,signal)

        return {'detector': self.detector_signal, 'snr': self.detector_sn, 'saturation': self.saturation, \
                'ngroups_map': self.ngroups_map}, self.bg_pix

    def _project(self, signal, extracted):

        return signal.grid, extracted['wavelength'], signal.pixgrid_list[0]


class SlitlessProjection(Projection):

    def _2d(self, extracted, saturation, signal):

        Projection._2d(self,extracted,saturation,signal)

        if self.signal.dispersion_axis == 'y':  # where projection is slitless and axis is y
            # need to rotate these 90 deg clockwise to match our normal axis orientation. np.rot90 only works CCW
            # so we need to rotate 3 times.
            det_sn = np.rot90(self.detector_sn, 3)
            det_signal = np.rot90(self.detector_signal, 3)
            det_sat = np.rot90(self.saturation, 3)
            det_groups = np.rot90(self.ngroups_map, 3)
        else:
            det_signal = self.detector_signal
            det_sn = self.detector_sn
            det_sat = self.saturation
            det_groups = self.ngroups_map

        return {'detector_unrotated': self.detector_signal, 'snr_unrotated': self.detector_sn, \
                'saturation_unrotated': self.saturation, 'ngroups_map_unrotated': self.ngroups_map, \
                'detector':det_signal, 'snr':det_sn, 'saturation':det_sat, 'ngroups_map': det_groups}, \
               self.bg_pix

    def _project(self, signal, extracted):
        '''
        This function handles (almost) all projection-type dependent factors (which are independent of strategy)
        This primarily includes wave_pix and pix_grid.

        The 2D projection types for slitless and multiorder are handled in _2d()
        '''
        # this is the spatial grid for the calculation. It needs to be added as an attribute for the
        # common/scene/coordinates tests to work.
        grid = signal.grid
        wave_pix = extracted['wavelength']
        if signal.dispersion_axis == 'y':
            wave_pix = wave_pix[::-1]
        # this is the detector plane pixel grid. for most modes we just grab and use it directly.
        # however, slitless we want to redefine it to be spatial on both axes.
        orig_grid = signal.pixgrid_list[0]
        # detector plane gets rotated depending on dispersion_axis so adjust accordingly
        if signal.dispersion_axis == 'x':
            pix_grid = Grid(grid.xsamp, orig_grid.ysamp, orig_grid.nx, orig_grid.ny)
        else:
            pix_grid = Grid(grid.ysamp, orig_grid.xsamp, orig_grid.ny, orig_grid.nx)

        return grid, wave_pix, pix_grid


class MultiorderProjection(SlitlessProjection):

    def _2d(self, extracted, saturation, signal):
        Projection._2d(self,extracted,saturation,signal)
        
        return {'detector_unrotated': self.detector_signal, 'snr_unrotated': self.detector_sn, \
                'saturation_unrotated': self.saturation, 'ngroups_map_unrotated': self.ngroups_map, \
                'detector':np.rot90(self.detector_signal), 'snr':np.rot90(self.detector_sn),
                'saturation':np.rot90(self.saturation), 'ngroups_map': np.rot90(self.ngroups_map)}, self.bg_pix

    def _project(self, signal, extracted):
        '''
        This function handles (almost) all projection-type dependent factors (which are independent of strategy)
        This primarily includes wave_pix and pix_grid.

        The 2D projection types for slitless and multiorder are handled in _2d()
        '''
        # this is the spatial grid for the calculation. It needs to be added as an attribute for the
        # common/scene/coordinates tests to work.
        grid = signal.grid
        # this is the wavelength sampling on the detector.
        wave_pix = extracted['wavelength']  # this is already a 1D np.array

        orig_grid = signal.pixgrid_list[0]
        pix_grid = IrregularGrid(np.arange(orig_grid.nx,0,-1), np.arange(orig_grid.ny))

        return grid, wave_pix, pix_grid
