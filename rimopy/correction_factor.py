"""Correction Factor

This module contains the coefficients of the RIMO correction factor (RCF).

This value corrects the simulated extra terrestrial lunar irradiance in order
to use it with photometers.

See "Roman et al., 2020: Correction of a lunar-irradiance model for aerosol optical depth
retrieval and comparison with a star photometer" for mor information.

It exports the following classes:
    * CorrectionParams - DataClass that contains the estimated coefficients of the
        RCF for a wavelength.

It exports the following functions:

    * get_correction_params - returns the RCF coefficients estimated for a wavelength
"""

from dataclasses import dataclass
from typing import List, Iterable

from numpy.typing import NDArray
import numpy as np

@dataclass
class CorrectionParams:
    """
    DataClass that contains the estimated coefficients of the RCF for a wavelength.

    Attributes
    ----------
    a_coeff : np.array of float
        RCF coefficient 'a'
    b_coeff : np.array of float
        RCF coefficient 'b'
    c_coeff : np.array of float
        RCF coefficient 'c'
    """
    a_coeff: NDArray[np.float32]
    b_coeff: NDArray[np.float32]
    c_coeff: NDArray[np.float32]

def _get_corrected_wavelengths() -> List[float]:
    """Gets all wavelengths (in nanometers) presented in the RCF model

    Returns
    -------
    list of float
        A list of floats that are the wavelengths in nanometers, in order
    """
    return [340, 380, 440, 500, 675, 870, 935, 1020, 1640]

def _get_all_correction_params() -> List[List[float]]:
    """Gets all RCF coefficients

    Returns
    -------
    list of list of float
        A list containing multiple list of floats. Each sublist is the list of coefficients
        for a wavelength
    """
    return [[1.186, -2.35 * 10**-2, 1.92 * 10**-1], [1.082, -4.17 * 10**-3, 7.10 * 10**-2],
            [1.062, -5.35 * 10**-4, 1.14 * 10**-2], [1.078, -8.93 * 10**-4, 1.11 * 10**-2],
            [1.092, -4.50 * 10**-4, 1.38 * 10**-2], [1.075, -2.05 * 10**-3, 1.37 * 10**-2],
            [1.071, -2.41 * 10**-3, 1.36 * 10**-2], [1.035, 5.55 * 10**-3, 2.79 * 10**-2],
            [1.047, -1.25 * 10**-3, 2.26 * 10**-2]]

def _get_all_as() -> List[float]:
    """Gets all 'a' RCF coefficients

    Returns
    -------
    list of float
        A list containing all 'a' coefficients in wavelength order
    """
    return list(map(lambda x: x[0], _get_all_correction_params()))

def _get_all_bs() -> List[float]:
    """Gets all 'b' RCF coefficients

    Returns
    -------
    list of float
        A list containing all 'b' coefficients in wavelength order
    """
    return list(map(lambda x: x[1], _get_all_correction_params()))

def _get_all_cs() -> List[float]:
    """Gets all 'c' RCF coefficients

    Returns
    -------
    list of float
        A list containing all 'c' coefficients in wavelength order
    """
    return list(map(lambda x: x[2], _get_all_correction_params()))

def _get_interpolated_correction_params(wavelengths_nm: Iterable[float]) -> 'CorrectionParams':
    """Estimate the RCF params with interpolation

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths (in nanometers) of which one wants to linearly interpolate the RCF params

    Returns
    -------
    'CorrectionParams'
        Estimated correction params
    """
    x_values = _get_corrected_wavelengths()
    wavelengths_nm = np.where(wavelengths_nm < x_values[0], x_values[0], wavelengths_nm)
    wavelengths_nm = np.where(wavelengths_nm > x_values[-1], x_values[-1], wavelengths_nm)
    all_as = _get_all_as()
    a_coeff = np.interp(wavelengths_nm, x_values, all_as)
    all_bs = _get_all_bs()
    b_coeff = np.interp(wavelengths_nm, x_values, all_bs)
    all_cs = _get_all_cs()
    c_coeff = np.interp(wavelengths_nm, x_values, all_cs)
    return CorrectionParams(a_coeff, b_coeff, c_coeff)

def get_correction_params(wavelengths_nm: Iterable[float]) -> 'CorrectionParams':
    """Gets the RCF correction parameters for a specific wavelength in nanometers

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths (in nanometers) of which one wants to estimate the RCF params

    Returns
    -------
    'CorrectionParams'
        Estimated correction params
    """
    return _get_interpolated_correction_params(wavelengths_nm)
