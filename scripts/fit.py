"""
Contains various functions for fitting Gaussian curves.
"""
import numpy as np


def gaussian(x, *params):
    """
    This function will generate a curve that is the sum of several Gaussians.
    The first value in the params array should be the y-axis offset
    The rest of the values in params should be in "sets" of 3,
    as each curve needs 3 guess parameters, A, x0, and c.
    """
    y0 = params[0]

    gSum = 0.0
    totalFits = int((len(params) - 1) / 3)

    for j in range(totalFits):
        A = params[(j * 3) + 1]
        x0 = params[(j * 3) + 2]
        c = params[(j * 3) + 3]

        gSum += y0 + A * np.exp(-((x - x0) ** 2) / (2 * c * c))

    return gSum


def super_gaussian(x, *params):
    """
    This function will return a super-Gaussian.
    The variable n can be modified to change the "squareness" of the curve.
    A value of n = 2 would be a normal Gaussian curve.
    """
    y0 = params[0]
    A = params[1]
    x0 = params[2]
    c = params[3]
    n = 6

    return y0 + A * np.exp(-(((x - x0) / (np.sqrt(2) * c)) ** n))


def fwhm(x0, c, n=2):
    """
    Calculates the two FWHM loci, and returns the distance between them.
    Recall that n = 2 for a Gaussian.
    """
    x1 = x0 + (np.log(2)) ** (1 / n) * (2 ** (1 / 2)) * c
    x2 = x0 - (np.log(2)) ** (1 / n) * (2 ** (1 / 2)) * c
    return abs(x1 - x2)
