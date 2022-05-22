"""
Dedicated code for reading the SIAF. Adding it here in case
pysiaf does not officially support Roman.
"""

from collections import OrderedDict
import sys

if sys.version_info < (3, 9):
    import importlib_resources
else:
    import importlib.resources as importlib_resources

import ast
import numpy as np
from xml.etree import ElementTree as ET

from pysiaf import aperture
from pysiaf.siaf import ApertureCollection


class RomanSiaf(ApertureCollection):
    """
    Class for reading and handling Roman SIAF information. This is based on
    pysiaf base classes, so all of the pysiaf functionality for, e.g., JWST
    SIAF objects is included here as well.
    """

    def __init__(self, siaf_file=None):

        if not siaf_file:
            with importlib_resources.path('soc_roman_tools.resources.data',
                                          'roman_siaf.xml') as siaf_file:
                self.siaf_file = siaf_file
        else:
            self.siaf_file = siaf_file

        self.apertures = self.read_roman_siaf()
        self.observatory = 'Roman'

    def read_roman_siaf(self):
        """
        Purpose
        -------
        Read the SIAF file.

        Inputs
        ------
        None

        Returns
        -------
        apertures (collections.OrderedDict)
            An ordered dictionary of pysiaf.aperture objects containing each
            aperture from the Roman SIAF file that was read.
        """

        apertures = OrderedDict()
        tree = ET.parse(self.siaf_file)
        for entry in tree.getroot().iter('SiafEntry'):
            roman_aperture = RomanAperture()
            apertype = None
            for node in entry:
                if node.tag == 'AperType':
                    apertype = node.text
                if (node.tag in aperture.ATTRIBUTES_THAT_CAN_BE_NONE) and \
                        (node.text is None):
                    value = node.text
                elif node.tag in aperture.INTEGER_ATTRIBUTES:
                    try:
                        value = int(node.text)
                    except (TypeError, ValueError) as e:
                        if node.tag == 'DetSciYAngle':
                            value = np.int(float(node.text))
                        elif not node.text:
                            value = None
                        else:
                            raise TypeError
                elif node.tag in aperture.STRING_ATTRIBUTES:
                    value = node.text
                elif node.tag in aperture.FLOAT_ATTRIBUTES:
                    value = float(node.text)
                # If it has children (which we can test by a simple boolean),
                # then we need to get things out of it. The only time this will
                # happen is the containers for the prism and grism parameters.
                # Each child will be a Position container.
                elif node.tag:
                    if apertype == 'FULLSCA':
                        data = list()
                        for child in node:
                            data.append(SpecPars(child))
                        value = data
                else:
                    try:
                        value = float(node.text)
                    except TypeError:
                        print('{}: {}'.format(node.tag, node.text))
                        raise TypeError
                setattr(roman_aperture, node.tag, value)

            apertures[roman_aperture.AperName] = roman_aperture

        return apertures


class RomanAperture(aperture.Aperture):
    """
    A class that defines the accepted aperture types for
    the Roman SIAF. This is built upon the pysiaf.aperture.Aperture
    base class.
    """

    _accepted_aperture_types = ['FULLSCA', 'COMPOUND']

    def __init__(self):
        super(RomanAperture, self).__init__()
        self.observatory = 'Roman'


class SpecPars:
    """
    Class to contain the aXe-like parameterization of spectroscopic
    parameters for wavelength-to-pixel conversion for Roman WFI
    spectroscopic guiding.

    Note that pixels are in the SCIENCE frame!
    """

    def __init__(self, spec_node):
        """
        Inputs
        ------
        spec_node (`xml.etree`):
            XML object containing the spectroscopic mode SIAF parameters.

        Returns
        -------
        None
        """

        # General info
        self.position = list()

        # X limits
        self.x_min = list()
        self.x_max = list()

        # Y limits
        self.y_min = list()
        self.y_max = list()

        # Spectroscopic parameters
        self.blue_x = list()
        self.blue_y = list()
        self.red_x = list()
        self.red_y = list()
        self.red_wave = list()
        self.blue_wave = list()

        for position in spec_node:
            if 'SpectralElement' in position.tag:
                self.mode = spec_node.attrib(['SpectralElement'])
            for child in position:
                if 'GridPosition' in child.tag:
                    self.position.append(int(child.tag))
                if 'XMin' in child.tag:
                    self.x_min.append(int(child.text))
                if 'XMax' in child.tag:
                    self.x_max.append(int(child.text))
                if 'YMin' in child.tag:
                    self.y_min.append(int(child.text))
                if 'YMax' in child.tag:
                    self.y_max.append(int(child.text))
                if 'BlueDeltaX' in child.tag:
                    self.blue_x.append(float(child.text))
                if 'BlueDeltaY' in child.tag:
                    self.blue_y.append(float(child.text))
                if 'RedDeltaX' in child.tag:
                    self.red_x.append(float(child.text))
                if 'RedDeltaY' in child.tag:
                    self.red_y.append(float(child.text))
                if 'BlueWave20' in child.tag:
                    self.blue_wave.append(float(child.text))
                if 'RedWave20' in child.tag:
                    self.blue_wave.append(float(child.text))
                else:
                    pass


def get_distortion_coeffs(aperture_name, siaf_file=None, inverse=False):
    """
    Purpose
    -------
    Parse the SIAF XML tree and pull out the distortion coefficients
    for a particular aperture name.

    Inputs
    ------
    aperture_name (string):
        Name of the aperture in the SIAF XML tree from which to pull the
        distortion coefficients.

    siaf_file (string; optional; default=None:
        Path to the SIAF file to be read. If None, read the default SIAF that
        comes with soc_roman_tools.

    inverse (boolean; optional; default=False):
        If True, return the inverse distortion coefficients.

    Returns
    -------
    x_coeffs (`~numpy.ndarray`):
        Array of polynomial coefficients describing the geometric distortion
        in the X direction.

    y_coeffs (`~numpy.ndarray`):
        Array of polynomial coefficients describing the geometric distortion
        in the Y direction.
    """

    if not siaf_file:
        with importlib_resources.path('soc_roman_tools.resources.data',
                                      'roman_siaf.xml') as sf:
            siaf_file = sf

    siaf_xml = ET.parse(siaf_file)
    siaf_root = siaf_xml.getroot()

    x_coeffs = {}
    y_coeffs = {}

    for child in siaf_root:
        if not inverse:
            if child[1].text == aperture_name:
                # X = indices 33 to 53
                # Y = indices 54 to 74
                for i in np.arange(33, 54, 1):
                    name = child[i].tag
                    j = ast.literal_eval(name[-2])
                    k = ast.literal_eval(name[-1])
                    value = float(child[i].text)
                    x_coeffs[f'c{j - k}_{k}'] = value
                for i in np.arange(54, 75, 1):
                    name = child[i].tag
                    j = ast.literal_eval(name[-2])
                    k = ast.literal_eval(name[-1])
                    value = float(child[i].text)
                    y_coeffs[f'c{j - k}_{k}'] = value
        else:
            if child[1].text == aperture_name:
                # X = indices 75 to 95
                # Y = indices 96 to 116
                for i in np.arange(75, 95, 1):
                    name = child[i].tag
                    j = ast.literal_eval(name[-2])
                    k = ast.literal_eval(name[-1])
                    value = float(child[i].text)
                    x_coeffs[f'c{j - k}_{k}'] = value
                for i in np.arange(96, 116, 1):
                    name = child[i].tag
                    j = ast.literal_eval(name[-2])
                    k = ast.literal_eval(name[-1])
                    value = float(child[i].text)
                    y_coeffs[f'c{j - k}_{k}'] = value

    return x_coeffs, y_coeffs
