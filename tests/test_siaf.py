import numpy as np

from soc_roman_tools.siaf import siaf


def test_default_siaf():
    # Just make sure this doesn't throw an exception.
    sfile = siaf.RomanSiaf().read_roman_siaf()


def test_coordinate_transforms():
    # Very basic coordinate transforms for now.
    # Using WFI01.
    sfile = siaf.RomanSiaf().read_roman_siaf()
    x1_det, y1_det = 1, 1
    x1_sci, y1_sci = -3, 4092

    det_sci_1 = sfile['WFI01_FULL'].det_to_sci(x1_det, y1_det)
    sci_det_1 = sfile['WFI01_FULL'].sci_to_det(x1_sci, y1_sci)

    # Move to the center of the detector to check the ideal frame.
    x2_det, y2_det = 2048.5, 2048.5
    x2_sci, y2_sci = 2044.5, 2044.5
    x2_idl, y2_idl = 0, 0

    det_sci_2 = sfile['WFI01_FULL'].det_to_sci(x2_det, y2_det)
    sci_det_2 = sfile['WFI01_FULL'].sci_to_det(x2_sci, y2_sci)

    det_idl_2 = sfile['WFI01_FULL'].det_to_idl(x2_det, y2_det)
    idl_det_2 = sfile['WFI01_FULL'].idl_to_det(x2_idl, y2_idl)

    sci_idl_2 = sfile['WFI01_FULL'].sci_to_idl(x2_sci, y2_sci)
    idl_sci_2 = sfile['WFI01_FULL'].idl_to_sci(x2_idl, y2_idl)

    # Now check all the transformations.
    np.testing.assert_allclose(det_sci_1[0], x1_sci, atol=0.01)
    np.testing.assert_allclose(det_sci_1[1], y1_sci, atol=0.01)

    np.testing.assert_allclose(sci_det_1[0], x1_det, atol=0.01)
    np.testing.assert_allclose(sci_det_1[1], y1_det, atol=0.01)

    np.testing.assert_allclose(det_sci_2[0], x2_sci, atol=0.01)
    np.testing.assert_allclose(det_sci_2[1], y2_sci, atol=0.01)

    np.testing.assert_allclose(sci_det_2[0], x2_det, atol=0.01)
    np.testing.assert_allclose(sci_det_2[1], y2_det, atol=0.01)

    np.testing.assert_allclose(det_idl_2[0], x2_idl, atol=0.01)
    np.testing.assert_allclose(det_idl_2[1], y2_idl, atol=0.01)

    np.testing.assert_allclose(idl_det_2[0], x2_det, atol=0.01)
    np.testing.assert_allclose(idl_det_2[1], y2_det, atol=0.01)

    np.testing.assert_allclose(sci_idl_2[0], x2_idl, atol=0.01)
    np.testing.assert_allclose(sci_idl_2[1], y2_idl, atol=0.01)

    np.testing.assert_allclose(idl_sci_2[0], x2_sci, atol=0.01)
    np.testing.assert_allclose(idl_sci_2[1], y2_sci, atol=0.01)
