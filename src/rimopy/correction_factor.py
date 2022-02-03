"""Correction Factor

This module contains the coefficients of the RIMO correction factor (RCF).

It exports the foollowing functions:

    * getCorrectionParams - returns the RCF coefficients estimated for a wavelength
"""

from argparse import ArgumentError
from typing import List
from scipy.interpolate import interp1d

class CorrectionParams:
    """
    DataModel that contains the estimated coefficients of the RCF for a wavelength.

    Attributes
    ----------
    a : float
        RCF coefficient 'a'
    b : float
        RCF coefficient 'b'
    c : float
        RCF coefficient 'c'
    """
    __slots__ = ['a', 'b', 'c']

    def __init__(self, a: float, b: float, c: float):
        """
        Parameters
        ----------
        a : float
            RCF coefficient 'a'
        b : float
            RCF coefficient 'b'
        c : float
            RCF coefficient 'c'
        """
        self.a = a
        self.b = b
        self.c = c

def _getCorrectedWavelengths() -> List[float]:
    """Gets all wavelengths (in nanometers) presented in the RCF model

    Returns
    -------
    list of float
        A list of floats that are the wavelengths in nanometers, in order
    """
    return [340, 380, 440, 500, 675, 870, 935, 1020, 1640]

def _getAllCorrectionParams() -> List[List[float]]:
    """Gets all RCF coefficients

    Returns
    ------- 
    list of list of float
        A list containing multiple list of floats. Each sublist is the list of coefficients for a wavelength
    """
    return [ [1.186, -2.35 * 10**-2, 1.92 * 10**-1 ], [1.082, -4.17 * 10**-3, 7.10 * 10**-2], [1.062, -5.35 * 10**-4, 1.14 * 10**-2],
        [1.078, -8.93 * 10**-4, 1.11 * 10**-2], [1.092, -4.50 * 10**-4, 1.38 * 10**-2], [1.075, -2.05 * 10**-3, 1.37 * 10**-2],
        [1.071, -2.41 * 10**-3, 1.36 * 10**-2], [1.035, 5.55 * 10**-3, 2.79 * 10**-2], [1.047, -1.25 * 10**-3, 2.26 * 10**-2] ]

def _getAllAs() -> List[float]:
    """Gets all 'a' RCF coefficients

    Returns
    ------- 
    list of float
        A list containing all 'a' coefficients in wavelength order
    """
    return list(map(lambda x : x[0], _getAllCorrectionParams()))

def _getAllBs() -> List[float]:
    """Gets all 'b' RCF coefficients

    Returns
    ------- 
    list of float
        A list containing all 'b' coefficients in wavelength order
    """
    return list(map(lambda x : x[1], _getAllCorrectionParams()))

def _getAllCs() -> List[float]:
    """Gets all 'c' RCF coefficients

    Returns
    ------- 
    list of float
        A list containing all 'c' coefficients in wavelength order
    """
    return list(map(lambda x : x[2], _getAllCorrectionParams()))

def _getInterpolatedCorrectionParams(wavelength_nm: float, kind='linear') -> 'CorrectionParams':
    """Estimate the RCF params with interpolation

    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which one wants to interpolate the RCF params
    kind: str
        Kind of interpolation performed. Specifies the kind of interpolation as a string or as an integer
        specifying the order of the spline interpolator to use. The string has to be one of 'linear', 'nearest',
        'nearest-up', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', or 'next'. Default is 'linear'.

    Raises
    ------
    ArgumentError
        If the argument kind does not fit the requirements specified in the section Parameters.

    Returns
    ------- 
    'CorrectionParams'
        Estimated correction params
    """
    x = _getCorrectedWavelengths()
    if wavelength_nm < x[0]:
        # Is this the best solution?
        return getCorrectionParams(x[0])
    if wavelength_nm > x[-1]:
        # Is this the best solution?
        return getCorrectionParams(x[-1])
    if kind not in ['linear', 'nearest', 'nearest-up', 'zero', 'slinear', 'quadratic',
        'cubic', 'previous', 'next'] and not isinstance(kind, int):
        raise ArgumentError("%s is unsupported, use a valid argument", kind)
    ya = _getAllAs()
    fa = interp1d(x, ya, kind)
    a = fa(wavelength_nm).item()
    yb = _getAllBs()
    fb = interp1d(x, yb, kind)
    b = fb(wavelength_nm).item()
    yc = _getAllCs()
    fc = interp1d(x, yc, kind)
    c = fc(wavelength_nm).item()
    return CorrectionParams(a, b, c)

def getCorrectionParams(wavelength_nm: float) -> 'CorrectionParams':
    """Gets the RCF correction parameters for a specific wavelength in nanometers

    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which one wants to estimate the RCF params

    Returns
    ------- 
    'CorrectionParams'
        Estimated correction params
    """
    wvs = _getCorrectedWavelengths()
    if wavelength_nm in wvs:
        index = wvs.index(wavelength_nm)
        corr_params = _getAllCorrectionParams()
        return CorrectionParams(corr_params[index][0], corr_params[index][1], corr_params[index][2])
    return _getInterpolatedCorrectionParams(wavelength_nm)