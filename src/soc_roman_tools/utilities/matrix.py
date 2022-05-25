"""
Adapted from an excerpt of Colin Cox's polynomial.py code written on 2013-04-29.
"""

import scipy


def dpdx(a, x, y, order=5):
    """
    Purpose
    -------
    Compute the differential of a 2D polynomial with respect to the X variable.
    Assumes that coefficients are ordered as described in the JWST and Roman
    science instrument aperture files (SIAFs).

    Inputs
    ------
    a (iterable of floats):
        An iterable (list or array) of float values giving the coefficients of
        a 2D polynomial as a function of X and Y.

    x (float or array of floats):
        The X coordinate(s) (with respect to the reference pixel position)
        at which to evaluate the partial derivative.

    y (float or array of floats):
        The Y coordinate(s) (with respect to the reference pixel position)
        at which to evaluate the partial derivative.

    order (integer; default = 5):
        The order of the 2D polynomial. Note that for Roman WFI, this should
        always be 5.

    Returns
    -------
    partial_x (float):
        The partial derivative with respect to X.
    """

    partial_x = 0.0
    k = 1  # index for coefficients
    for i in range(1, order + 1):
        for j in range(i + 1):
            if i - j > 0:
                partial_x += (i - j) * a[k] * x**(i - j - 1) * y**j
            k += 1
    return partial_x


def dpdy(a, x, y, order=5):
    """
    Purpose
    -------
    Compute the differential of a 2D polynomial with respect to the Y variable.
    Assumes that coefficients are ordered as described in the JWST and Roman
    science instrument aperture files (SIAFs).

    Inputs
    ------
    b (iterable of floats):
        An iterable (list or array) of float values giving the coefficients of
        a 2D polynomial as a function of X and Y.

    x (float or array of floats):
        The X coordinate(s) (with respect to the reference pixel position)
        at which to evaluate the partial derivative.

    y (float or array of floats):
        The Y coordinate(s) (with respect to the reference pixel position)
        at which to evaluate the partial derivative.

    order (integer; default = 5):
        The order of the 2D polynomial. Note that for Roman WFI, this should
        always be 5.

    Returns
    -------
    partial_x (float):
        The partial derivative with respect to X.
    """

    partial_y = 0.0
    k = 1  # index for coefficients
    for i in range(1, order + 1):
        for j in range(i + 1):
            if j > 0:
                partial_y += j * a[k] * x**(i - j) * y**(j - 1)
            k += 1
    return partial_y


def jacob(a, b, x, y, order=5):
    """
    Purpose
    -------
    Compute the Jacobian determinant of a 2D polynomial. In principal,
    this is used to compute the area of each pixel on the sky given
    a polynomial that describes the geometric distortion.

    Note that the functions called (dpdx and dpdy) assume that the
    order of the polynomial coefficients is the order used in the
    JWST and Roman science instrument aperture files (SIAFs).

    Inputs
    ------
    a (iterable of floats):
        An iterable (list or array) of float values giving the coefficients of
        a 2D polynomial as a function of X and Y that fit the X pixel positions,
        i.e., X - Xsci = f(X, Y).

    b (iterable of floats):
        An iterable (list or array) of float values giving the coefficients of
        a 2D polynomial as a function of X and Y that fit the Y pixel positions,
        i.e., Y - Ysci = f(X, Y).

    x (float or iterable of floats):
        The X coordinate(s) (with respect to the reference pixel position)
        at which to evaluate the partial derivative.

    y (float or iterable of floats):
        The Y coordinate(s) (with respect to the reference pixel position)
        at which to evaluate the partial derivative.

    order (integer; default = 5):
        The order of the 2D polynomial. Note that for Roman WFI, this should
        always be 5.

    Returns
    -------
    jacobian (float or iterable of floats):
        The Jacobian determinant of the 2D polynomial evaluated at one or more
        input positions. The shape of the result will match the shape of the
        input x and y variables.

        This is the area on the sky of a pixel at position (x, y) given a
        geometric distortion described by a 2D polynomial with coefficients
        a and b as further described in the JWST and Roman SIAF documentation.
    """
    jacobian = dpdx(a, x, y, order=order) * dpdy(b, x, y, order=order) - \
               dpdx(b, x, y, order=order) * dpdy(a, x, y, order=order)
    jacobian = scipy.fabs(jacobian)
    return jacobian
