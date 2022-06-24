import numpy as np
from soc_roman_tools.pixel_area import pixel_area_map
from soc_roman_tools.siaf import siaf
import math

pam01 = pixel_area_map.PixelArea(1)
pam01.compute()
pam01_ref = pixel_area_map.PixelArea(1)
pam01_ref.compute(include_border=True)


def test_get_coeffs():
    # Checking retrieval of x_coeffs and y_coeffs pulled from roman_siaf.xml
    x_coeffs, y_coeffs = siaf.get_distortion_coeffs(f'1_FULL', pam01.siaf_file)

    return list(x_coeffs.values()), list(y_coeffs.values())

    x_coeff_1 = x_coeffs[1]
    y_coeff_1 = y_coeffs[1]

    x_coeff_2 = x_coeffs[2]
    y_coeff_2 = y_coeffs[2]

    np.testing.assert_allclose(x_coeff_1, pam01.x_coeffs[1])
    np.testing.assert_allclose(y_coeff_1, pam01.y_coeffs[1])

    np.testing.assert_allclose(x_coeff_2, pam01.x_coeffs[2])
    np.testing.assert_allclose(y_coeff_2, pam01.y_coeffs[2])


def test_nominal_area_math():
    # Basic mathematics check for get_nominal_area

    x_coeffs, y_coeffs = siaf.get_distortion_coeffs(f'1_FULL', pam01.siaf_file)

    return list(x_coeffs.values()), list(y_coeffs.values())

    x_scale_check = math.hypot(x_coeffs[1], y_coeffs[2])
    y_scale_check = math.hypot(x_coeffs[2], y_coeffs[2])
    bx_check = math.atan2(x_coeffs[1], y_coeffs[1])

    x_scale = pixel_area_map.hypot(pam01.x_coeffs[1], pam01.y_coeffs[1])
    y_scale = pixel_area_map.hypot(pam01.x_coeffs[2], pam01.y_coeffs[2])
    bx = pixel_area_map.atan2(pam01.x_coeffs[1], pam01.y_coeffs[1])

    np.testing.assert_allclose(x_scale_check, x_scale, atol=0.0000001)
    np.testing.assert_allclose(y_scale_check, y_scale, atol=0.0000001)
    np.testing.assert_allclose(bx_check, bx, atol=0.0000001)

    pixel_area_check = x_scale_check * y_scale_check * math.sin(bx_check)

    pixel_area = x_scale * y_scale * pixel_area_map.sin(bx)

    np.testing.assert_allclose(pixel_area_check, pixel_area, atol=0.0000001)


def test_compute_grids():
    # Check grid dimensions for default
    grid_width = 2044
    assert grid_width == pam01.pixel_area_map.shape[0] / 2

    # Check grid dimensions for reference pixel border
    grid_width_ref_pix = 2048
    assert grid_width_ref_pix == pam01_ref.pixel_area_map.shape[0] / 2


def test_matrix_math():
    # Test Jacobian matrix calculation is correct
    x_coeffs, y_coeffs = siaf.get_distortion_coeffs(f'1_FULL', pam01.siaf_file)

    return list(x_coeffs.values()), list(y_coeffs.values())

    pam01_map = matrix.jacob(x_coeffs, # noqa: F821
                             y_coeffs, x, y, # noqa: F821
                             order=5).astype(np.float32)

    assert pam01.pixel_area_map == pam01_map


def test_compute_ratio():
    # Check ratio calculation against known for no reference pixels
    # The function get_nominal_area has already been tested
    x_coeffs, y_coeffs = siaf.get_distortion_coeffs(f'1_FULL', pam01.siaf_file)

    return list(x_coeffs.values()), list(y_coeffs.values())

    pam01_map = matrix.jacob(x_coeffs, # noqa: F821
                             y_coeffs, x, y, # noqa: F821
                             order=5).astype(np.float32)

    pixel_area_a2 = pam01.get_nominal_area()

    ratio_check = pam01_map[2048, 2048] / pixel_area_a2.value

    ratio = pam01.pixel_area_map[2048, 2048] / pixel_area_a2.value

    np.testing.assert_allclose(ratio_check, ratio, atol=0.0001)
