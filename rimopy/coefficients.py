"""Coefficients

This module contains the coefficient and wavelength data from the ROLO model.
Currently it only contains the first 9 wavelengths.

It exports the following functions:

    * get_wavelengths - returns all wavelengths present in the ROLO model
    * get_all_coefficients_a - returns all 'a' coefficients for all wavelengths
    * get_all_coefficients_b - returns all 'b' coefficients for all wavelengths
    * get_all_coefficients_d - returns all 'd' coefficients for all wavelengths
    * get_coefficients_a - returns the 'a' coefficients for a specific wavelength index
    * get_coefficients_b - returns the 'b' coefficients for a specific wavelength index
    * get_coefficients_d - returns the 'd' coefficients for a specific wavelength index
    * get_coefficients_c - returns the 'c' coefficients
    * get_coefficients_p - returns the 'p' coefficients
    * get_apollo_coefficients - returns all Apollo adjusting coefficients
"""
from dataclasses import dataclass
from typing import List, Dict, Iterable
import csv
import pkgutil
from io import StringIO

import numpy as np
from numpy.typing import NDArray

@dataclass
class _CoefficientsWln():
    """
    Coefficients data for a wavelength. It includes only the a, b and d coefficients.

    Attributes
    ----------
    a_coeffs : tuple of 4 floats, corresponding to coefficients a0, a1, a2, and a3
    b_coeffs : tuple of 3 floats, corresponding to coefficients b1, b2, and b3
    d_coeffs : tuple of 3floats, corresponding to coefficients d1, d2, and d3
    """
    __slots__ = ['a_coeffs', 'b_coeffs', 'd_coeffs']

    def __init__(self, coeffs: List[float]):
        """
        Parameters
        ----------
        coeffs : list of float
            List of floats consisting of all coefficients. In order: a0, a1, a2, a3, b1, b2, b3,
            d1, d2 and d3.

        Returns
        -------
        _CoefficientsWln
            Instance of _Coefficients with the correct data
        """
        self.a_coeffs = (coeffs[0], coeffs[1], coeffs[2], coeffs[3])
        self.b_coeffs = (coeffs[4], coeffs[5], coeffs[6])
        self.d_coeffs = (coeffs[7], coeffs[8], coeffs[9])

def _get_coefficients_data() -> Dict[float, '_CoefficientsWln']:
    """Returns all variable coefficients (a, b and d) for all wavelengths

    Returns
    -------
    A dict that has the wavelengths as keys (float), and as values the _CoefficientsWln associated
    to the wavelength.
    """
    coeff_bytes = pkgutil.get_data(__name__, 'data/coefficients.csv')
    coeff_string = coeff_bytes.decode()
    file = StringIO(coeff_string)
    csvreader = csv.reader(file)
    next(csvreader) # Discard the header
    data = {}
    for row in csvreader:
        coeffs = []
        for i in range(1, 11):
            coeffs.append(float(row[i]))
        data[float(row[0])] = _CoefficientsWln(coeffs)
    file.close()
    return data

def get_wavelengths() -> List[float]:
    """Gets all wavelengths present in the model, in nanometers

    Returns
    -------
    list of float
        A list of floats that are the wavelengths in nanometers, in order
    """
    coeffs = _get_coefficients_data()
    return list(coeffs.keys())

def get_all_coefficients_a() -> List[List[float]]:
    """Gets all 'a' coefficients

    Returns
    -------
    list of list of float
        A list containing multiple list of floats. Each sublist is the list of 'a' coefficients
        for a wavelength
    """
    coeffs = _get_coefficients_data()
    return [elem.a_coeffs for elem in coeffs.values()]

def get_all_coefficients_b() -> List[List[float]]:
    """Gets all 'b' coefficients

    Returns
    -------
    list of list of float
        A list containing multiple list of floats. Each sublist is the list of 'b' coefficients
        for a wavelength
    """
    coeffs = _get_coefficients_data()
    return [elem.b_coeffs for elem in coeffs.values()]

def get_all_coefficients_d() -> List[List[float]]:
    """Gets all 'd' coefficients

    Returns
    -------
    list of list of float
        A list containing multiple list of floats. Each sublist is the list of 'd' coefficients
        for a wavelength
    """
    coeffs = _get_coefficients_data()
    return [elem.d_coeffs for elem in coeffs.values()]

def get_coefficients_a(wavelengths_nm: NDArray[np.float32]) -> NDArray[np.float32]:
    """Gets all 'a' coefficients for a concrete wavelength

    Parameters
    ----------
    wavelength_nm : array of float
        Wavelengths in nanometers from which one wants to obtain the coefficients.

    Returns
    -------
    np.array of float
        An array containing the 'a' coefficients for the wavelengths.
        It contains 4 arrays, one per coeff, and each inner array contains one value per given wavelength.
    """
    a_coeffs = get_all_coefficients_a()
    wvs = get_wavelengths()
    wavelengths_nm = np.where(wavelengths_nm < wvs[0], wvs[0], wavelengths_nm)
    wavelengths_nm = np.where(wavelengths_nm > wvs[-1], wvs[-1], wavelengths_nm)
    a_0 = np.interp(wavelengths_nm, wvs, [elem[0] for elem in a_coeffs])
    a_1 = np.interp(wavelengths_nm, wvs, [elem[1] for elem in a_coeffs])
    a_2 = np.interp(wavelengths_nm, wvs, [elem[2] for elem in a_coeffs])
    a_3 = np.interp(wavelengths_nm, wvs, [elem[3] for elem in a_coeffs])
    a_coeffs_concrete = np.array([a_0, a_1, a_2, a_3])
    return a_coeffs_concrete

def get_coefficients_b(wavelengths_nm: NDArray[np.float32]) -> NDArray[np.float32]:
    """Gets all 'b' coefficients for a concrete wavelength

    Parameters
    ----------
    wavelengths_nm : array of float
        Wavelength in nanometers from which one wants to obtain the coefficients.

    Returns
    -------
    np.array of float
        An array containing the 'b' coefficients for the wavelengths.
        It contains 3 arrays, one per coeff, and each inner array contains one value per given wavelength.
    """
    b_coeffs = get_all_coefficients_b()
    wvs = get_wavelengths()
    wavelengths_nm = np.where(wavelengths_nm < wvs[0], wvs[0], wavelengths_nm)
    wavelengths_nm = np.where(wavelengths_nm > wvs[-1], wvs[-1], wavelengths_nm)
    b_0 = np.interp(wavelengths_nm, wvs, [elem[0] for elem in b_coeffs])
    b_1 = np.interp(wavelengths_nm, wvs, [elem[1] for elem in b_coeffs])
    b_2 = np.interp(wavelengths_nm, wvs, [elem[2] for elem in b_coeffs])
    b_coeffs_concrete = np.array([b_0, b_1, b_2])
    return b_coeffs_concrete

def get_coefficients_d(wavelengths_nm: NDArray[np.float32]) -> NDArray[np.float32]:
    """Gets all 'd' coefficients for a concrete wavelength

    Parameters
    ----------
    wavelengths_nm : array of float
        Wavelengths in nanometers from which one wants to obtain the coefficients.

    Returns
    -------
    np.array of float
        An array containing the 'd' coefficients for the wavelengths.
        It contains 3 arrays, and each inner array contains one value per given wavelength
    """
    d_coeffs = get_all_coefficients_d()
    wvs = get_wavelengths()
    wavelengths_nm = np.where(wavelengths_nm < wvs[0], wvs[0], wavelengths_nm)
    wavelengths_nm = np.where(wavelengths_nm > wvs[-1], wvs[-1], wavelengths_nm)
    d_0 = np.interp(wavelengths_nm, wvs, [elem[0] for elem in d_coeffs])
    d_1 = np.interp(wavelengths_nm, wvs, [elem[1] for elem in d_coeffs])
    d_2 = np.interp(wavelengths_nm, wvs, [elem[2] for elem in d_coeffs])
    d_coeffs_concrete = np.array([d_0, d_1, d_2])
    return d_coeffs_concrete

def get_coefficients_c() -> NDArray[np.float32]:
    """Gets all 'c' coefficients

    Returns
    -------
    np.array of float
        An array containing all 'c' coefficients
    """
    return np.array([0.00034115, -0.0013425, 0.00095906, 0.00066229])

def get_coefficients_p() -> NDArray[np.float32]:
    """Gets all 'p' coefficients

    Returns
    -------
    np.array of float
        An array containing all 'p' coefficients
    """
    return np.array([4.06054, 12.8802, -30.5858, 16.7498])

def get_apollo_coefficients() -> NDArray[np.float32]:
    """Coefficients used for the adjustment of the ROLO model using Apollo spectra.

    Returns
    -------
    np.array of float
        An array containing all Apollo coefficients
    """
    return np.array([1.0301, 1.0970, 0.9325, 0.9466, 1.0225, 1.0157, 1.0470, 1.0084, 1.0100, 1.0148,
            0.9843, 1.0134, 0.9329, 0.9849, 0.9994, 0.9957, 1.0059, 0.9618, 0.9561, 0.9796, 0.9568,
            0.9873, 1.0575, 1.0108, 0.9743, 1.0386, 1.0338, 1.0577, 1.0650, 1.0815, 0.8945, 0.9689])
