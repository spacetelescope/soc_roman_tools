"""
Tools for creating pixel area maps of the Roman WFI detectors.

Some of this code has been adapted from the JWST NIRCam code for
making their calibration reference files.
"""

import datetime
import getpass
import logging
from math import hypot, atan2, sin
import sys

if sys.version_info < (3, 9):
    import importlib_resources
else:
    import importlib.resources as importlib_resources

import asdf
from astropy import units as u, __version__ as astropy_version
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np

from roman_datamodels import stnode as rds, __version__ as rdm_version

from ..utilities import logging_functions
from ..utilities import matrix
from soc_roman_tools.siaf import siaf

# Automatic versioning
from soc_roman_tools.version import version as __version__


class PixelArea:
    """
    Class for computing a pixel area map for the Roman WFI detectors.

    Inputs
    ------
    detector (integer):
        WFI detector ID number.

    siaf_file (string; default=None):
        Path to the science instrument aperture file (SIAF) containing
        the geometric distortion polynomial coefficients. If None, then
        the code will look up the copy of the SIAF included in the
        package data.

    verbose (boolean; default=False):
        Optional argument to enable logging messages with additional
        information.

    Examples
    --------
    To compute the pixel area map of the science pixels (4088 x 4088) of
    detector WFI07 that does not include the reference pixel border:

    >>> from soc_roman_tools.pixel_area import pixel_area_map
    >>> pam07 = pixel_area_map.PixelArea(7)
    >>> pam07.compute()
    >>> pam07.save_asdf()
    """

    def __init__(self, detector, siaf_file=None, verbose=False):

        # Set up verbosity
        self.verbose = verbose
        if self.verbose:
            logging_functions.configure_logging(level=logging.INFO)
        else:
            logging_functions.configure_logging(level=logging.WARNING)

        # Some record-keeping attributes.
        if not siaf_file:
            with importlib_resources.path('soc_roman_tools.resources.data',
                                          'roman_siaf.xml') as sf:
                self.siaf_file = sf
        else:
            self.siaf_file = siaf_file
        self.detector = f'WFI{detector:02d}'

        # Get the distortion coefficients from the SIAF.
        self.x_coeffs, self.y_coeffs = self.get_coeffs()

        self.pixel_area_map = None
        self.nominal_pixel_area = None

        logging.info(f'Set up to create pixel area map for {self.detector}.')
        logging.info(f'soc_roman_tools version = {__version__}')
        logging.info(f'roman_datamodels version = {rdm_version}')
        logging.info(f'astropy version = {astropy_version}')
        logging.info(f'numpy version = {np.__version__}')

    @logging_functions.timeit
    def compute(self, include_border=False, refpix_area=False):
        """
        Purpose
        -------
        Perform the steps necessary to construct the pixel area map.

        Inputs
        ------
        include_border (boolean; default=False):
            Include the 4 pixel reference pixel border around the science
            pixel array. This makes the pixel area map have dimensions
            of 4096 x 4096. It is recommended to leave this set to False.

        refpix_area (boolean; default=False):
            If the reference pixel border is included, then this parameter
            indicates whether or not to set the area of the reference
            pixels to zero. A value of True will compute the area of the
            reference pixels. Default is False.

        Returns
        -------
        None
        """

        pixel_area_a2 = self.get_nominal_area()

        # Make a grid of pixel positions from from -det_size to +det_size.
        # Really this is half of the detector size, so that it's the full
        # detector, but centered at 0 (the reference pixel position).
        #
        # If for some reason the user wants to include the border of reference
        # pixels, we can do that here...might be useful if someone skipped
        # trimming them.
        if include_border:
            logging.info('Including reference pixel border in pixel area '
                         'map. Pixel area map will have dimensions 4096 x '
                         '4096.')
            logging.warning('!!')
            logging.warning('Calibrated Roman WFI data have the reference '
                            'pixel border trimmed, and dimensions of 4088 x '
                            '4088. Do not use this pixel area map unless '
                            'the calibrated data include the reference '
                            'pixel border.')
            logging.warning('!!')
            det_size = 2048
        else:
            det_size = 2044
        pixels = np.mgrid[-det_size:det_size, -det_size:det_size]
        y = pixels[0, :, :]
        x = pixels[1, :, :]

        self.pixel_area_map = matrix.jacob(self.x_coeffs,
                                           self.y_coeffs, x, y,
                                           order=5).astype(np.float32)

        # Sanity check.
        ratio = self.pixel_area_map[det_size, det_size] / pixel_area_a2.value
        logging.info(f'Jacobian determinant at (x, y) = (2044, 2044) is '
                     f'{self.pixel_area_map[det_size, det_size]} arcsec^2')
        logging.info(f'Nominal pixel area is {pixel_area_a2.value} '
                     f'arcsec^2')
        logging.info(f'Ratio (Jacobian / nominal) = {ratio}')

        # Normalize the pixel area map to the nominal pixel area.
        # Both are in units of arcseconds before the normalization.
        self.pixel_area_map /= pixel_area_a2.value

        # If the reference pixel border was included, check if we're
        # setting the area of the reference pixels to zero. This is the
        # default behavior...someone might override it.
        if include_border:
            if not refpix_area:
                logging.info(f'Reference pixel border was included, but '
                             f'refpix_area = {refpix_area}. Setting area '
                             f'of reference pixels to zero.')
                self.pixel_area_map[:4, :] = 0
                self.pixel_area_map[-4:, :] = 0
                self.pixel_area_map[:, :4] = 0
                self.pixel_area_map[:, -4:] = 0

        # Save the nominal pixel area in units of steradians.
        self.nominal_pixel_area = pixel_area_a2.to(u.sr)

    def save_asdf(self, filename=None, meta_data_override={}):
        """
        Purpose
        -------
        Write the pixel area map to an ASDF file using the
        AREA reference file datamodel.

        Inputs
        ------
        filename (string; default=None):
            Name of the output ASDF file. If None, then
            construct a file name of the form
            'roman_{detector}_YYYYMMDD_hhmmss_area.asdf'.
            For example:
                roman_wfi16_20220211_222140_area.asdf
            is a pixel area map of detector WFI16 that
            was made on February 11, 2022 at 22:21:40.

        meta_data_override (dictionary; default=None):
            A dictionary of values to override default
            entries in the output ASDF file meta data.

        Returns
        -------
        None

        Examples
        --------
        To generate a pixel area map of WFI16 and save the output
        to an ASDF file, while overriding the meta data to set the
        origin to STScI and the useafter to 2022-01-01 00:00:00:

        >>> from datetime import datetime
        >>> from soc_roman_tools.pixel_area import pixel_area_map
        >>> from astropy.time import Time
        >>> pam16 = pixel_area_map.PixelArea(16)
        >>> pam16.compute()
        >>> useafter = Time(datetime(2022, 1, 1, 0, 0, 0))
        >>> pam16.save_asdf(meta_data_override={'origin': 'STScI',
        >>>                                     'useafter': useafter})
        """

        if not filename:
            date = datetime.datetime.today().strftime('%Y%m%d_%H%M%S')
            filename = f'roman_{self.detector.lower()}_{date}_area.asdf'

        dm = rds.PixelareaRef()

        meta = {'reftype': 'AREA',
                'description': 'Roman WFI pixel area map.',
                'pedigree': 'GROUND',
                'telescope': 'ROMAN',
                'origin': getpass.getuser(),
                'author': f'wfitools version {__version__}',
                'useafter': Time(datetime.datetime(2020, 1, 1, 0, 0, 0)),
                'photometry':
                    {'pixelarea_arcsecsq':
                         self.nominal_pixel_area.to(u.arcsec * u.arcsec),
                     'pixelarea_steradians': self.nominal_pixel_area},
                'instrument':
                    {'optical_element': 'F158',
                     'detector': self.detector.upper(),
                     'name': 'WFI'}
                }

        # If the user wants to override any of these defaults, then do so now.
        meta.update(meta_data_override)

        logging.info('Saving ASDF file.')
        logging.info(f'ASDF meta data: {meta}')

        dm['data'] = self.pixel_area_map
        dm['meta'] = meta

        asdf_file = asdf.AsdfFile()
        asdf_file.tree = {'roman': dm}
        asdf_file.write_to(filename)

        logging.info(f'Pixel area map saved to file {filename}')

    def show_map(self, filename=None):
        """
        Purpose
        -------
        Display the pixel area map.

        Inputs
        ------
        filename (string; default=None):
            Name of the file to save the pixel area map image. If None,
            the image will be displayed on the screen instead.

        Returns
        -------
        None
        """

        fig = plt.figure('Pixel Area Map')
        ax = fig.add_subplot()
        img = ax.imshow(self.pixel_area_map, origin='lower',
                        vmin=np.min(self.pixel_area_map[self.pixel_area_map > 0]))
        ax.set_xlabel('X science coordinate (pixels)')
        ax.set_ylabel('Y science coordinate (pixels)')
        ax.set_title(f'Pixel Area Map for {self.detector}')
        plt.colorbar(img)

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    def get_nominal_area(self):
        """
        Purpose
        -------
        Compute the nominal pixel area at the reference position.

        Inputs
        ------
        None

        Returns
        -------
        pixel_area (`~astropy.units.Quantity`):
            The area of the nominal reference pixel in units of square
            arcseconds.
        """

        x_scale = hypot(self.x_coeffs[1], self.y_coeffs[1])
        y_scale = hypot(self.x_coeffs[2], self.y_coeffs[2])
        bx = atan2(self.x_coeffs[1], self.y_coeffs[1])

        pixel_area = x_scale * y_scale * sin(bx) * u.arcsec * u.arcsec

        return pixel_area

    def get_coeffs(self):
        """
        Purpose
        -------
        Get the geometric distortion polynomial coefficients from the SIAF.

        Inputs
        ------
        None

        Returns
        x_coeffs (list):
            Array of polynomial coefficients describing the geometric distortion
            in the X direction.

        y_coeffs (list):
            Array of polynomial coefficients describing the geometric distortion
            in the Y direction.
        """

        det_name = f'{self.detector}_FULL'
        x_coeffs, y_coeffs = siaf.get_distortion_coeffs(det_name, self.siaf_file)

        return list(x_coeffs.values()), list(y_coeffs.values())
