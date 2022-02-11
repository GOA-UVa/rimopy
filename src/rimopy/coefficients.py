"""Coefficients

This module contains the coefficient and wavelength data from the ROLO model.
Currently it only contains the first 9 wavelengths.

It exports the foollowing functions:

    * getWavelengths - returns all wavelengths present in the ROLO model
    * getAllCoefficientsA - returns all 'a' coefficients for all wavelengths
    * getAllCoefficientsB - returns all 'b' coefficients for all wavelengths
    * getAllCoefficientsD - returns all 'd' coefficients for all wavelengths
    * getCoefficientsA - returns the 'a' coefficients for a specific wavelength index
    * getCoefficientsB - returns the 'b' coefficients for a specific wavelength index
    * getCoefficientsD - returns the 'd' coefficients for a specific wavelength index
    * getCoefficientsC - returns the 'c' coefficients
    * getCoefficientsP - returns the 'p' coefficients
    * getApolloCoefficients - returns all Apollo adjusting coefficients
"""
from typing import List, Dict
from scipy.interpolate import interp1d
import csv
import pkgutil
from io import StringIO

_INTERPOLATION_TYPE = 'linear'

class _CoefficientsWln():
    """
    Coefficients data for a wavelength. It includes only the a, b and d coefficients.

    Attributes
    ----------
    a : tuple of 4 floats, corresponding to coefficients a0, a1, a2, and a3
    b : tuple of 3 floats, corresponding to coefficients b1, b2, and b3
    d : tuple of 3floats, corresponding to coefficients d1, d2, and d3
    """
    __slots__ = ['a', 'b', 'd']

    def __init__(self, cf: List[float]):
        """
        Parameters
        ----------
        cf : list of float
            List of floats consisting of all coefficients. In order: a0, a1, a2, a3, b1, b2, b3, d1, d2 and d3.
        """
        self.a = (cf[0], cf[1], cf[2], cf[3])
        self.b = (cf[4], cf[5], cf[6])
        self.d = (cf[7], cf[8], cf[9])
        

def _getCoefficientsData() -> Dict[float, '_CoefficientsWln']:
    """Returns all variable coefficients (a, b and d) for all wavelengths

    Returns
    -------
    A dict that has the wavelengths as keys (float), and as values the _CoefficientsWln associated to the wavelength.
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

def getWavelengths() -> List[float]:
    """Gets all wavelengths present in the model, in nanometers

    Returns
    -------
    list of float
        A list of floats that are the wavelengths in nanometers, in order
    """
    coeffs = _getCoefficientsData()
    return list(coeffs.keys())

def getAllCoefficientsA() -> List[List[float]]:
    """Gets all 'a' coefficients

    Returns
    ------- 
    list of list of float
        A list containing multiple list of floats. Each sublist is the list of 'a' coefficients for a wavelength
    """
    coeffs = _getCoefficientsData()
    return [elem.a for elem in coeffs.values()]

def getAllCoefficientsB() -> List[List[float]]:
    """Gets all 'b' coefficients

    Returns
    ------- 
    list of list of float
        A list containing multiple list of floats. Each sublist is the list of 'b' coefficients for a wavelength
    """
    coeffs = _getCoefficientsData()
    return [elem.b for elem in coeffs.values()]

def getAllCoefficientsD() -> List[List[float]]:
    """Gets all 'd' coefficients

    Returns
    ------- 
    list of list of float
        A list containing multiple list of floats. Each sublist is the list of 'd' coefficients for a wavelength
    """
    coeffs = _getCoefficientsData()
    return [elem.d for elem in coeffs.values()]

def getCoefficientsA(wavelength_nm: float) -> List[float]:
    """Gets all 'a' coefficients for a concrete wavelength

    Parameters
    ----------
    wavelength_nm : float
        Wavelength in nanometers from which one wants to obtain the coefficients.

    Returns
    ------- 
    list of float
        A list containing the 'a' coefficients for the wavelength
    """
    y = getAllCoefficientsA()
    wvs = getWavelengths()
    if wavelength_nm in wvs:
        return y[wvs.index(wavelength_nm)]
    if wavelength_nm < wvs[0]:
        return y[0]
    if wavelength_nm > wvs[-1]:
        return y[-1]
    f0 = interp1d(wvs, [elem[0] for elem in y], _INTERPOLATION_TYPE)
    f1 = interp1d(wvs, [elem[1] for elem in y], _INTERPOLATION_TYPE)
    f2 = interp1d(wvs, [elem[2] for elem in y], _INTERPOLATION_TYPE)
    f3 = interp1d(wvs, [elem[3] for elem in y], _INTERPOLATION_TYPE)
    y2 = [f0(wavelength_nm).item(), f1(wavelength_nm).item(), f2(wavelength_nm).item(), f3(wavelength_nm).item()]
    return y2

def getCoefficientsB(wavelength_nm: float) -> List[float]:
    """Gets all 'b' coefficients for a concrete wavelength

    Parameters
    ----------
    wavelength_nm : float
        Wavelength in nanometers from which one wants to obtain the coefficients.

    Returns
    ------- 
    list of float
        A list containing the 'b' coefficients for the wavelength
    """
    y = getAllCoefficientsB()
    wvs = getWavelengths()
    if wavelength_nm in wvs:
        return y[wvs.index(wavelength_nm)]
    if wavelength_nm < wvs[0]:
        return y[0]
    if wavelength_nm > wvs[-1]:
        return y[-1]
    f0 = interp1d(wvs, [elem[0] for elem in y], _INTERPOLATION_TYPE)
    f1 = interp1d(wvs, [elem[1] for elem in y], _INTERPOLATION_TYPE)
    f2 = interp1d(wvs, [elem[2] for elem in y], _INTERPOLATION_TYPE)
    y2 = [f0(wavelength_nm).item(), f1(wavelength_nm).item(), f2(wavelength_nm).item()]
    return y2

def getCoefficientsD(wavelength_nm: float) -> List[float]:
    """Gets all 'd' coefficients for a concrete wavelength

    Parameters
    ----------
    wavelength_nm : float
        Wavelength in nanometers from which one wants to obtain the coefficients.

    Returns
    ------- 
    list of float
        A list containing the 'd' coefficients for the wavelength
    """
    y = getAllCoefficientsD()
    wvs = getWavelengths()
    if wavelength_nm in wvs:
        return y[wvs.index(wavelength_nm)]
    if wavelength_nm < wvs[0]:
        return y[0]
    if wavelength_nm > wvs[-1]:
        return y[-1]
    f0 = interp1d(wvs, [elem[0] for elem in y], _INTERPOLATION_TYPE)
    f1 = interp1d(wvs, [elem[1] for elem in y], _INTERPOLATION_TYPE)
    f2 = interp1d(wvs, [elem[2] for elem in y], _INTERPOLATION_TYPE)
    y2 = [f0(wavelength_nm).item(), f1(wavelength_nm).item(), f2(wavelength_nm).item()]
    return y2

def getCoefficientsC() -> List[float]:
    """Gets all 'c' coefficients

    Returns
    ------- 
    list of float
        A list containing all 'c' coefficients
    """
    return [0.00034115, -0.0013425, 0.00095906, 0.00066229]

def getCoefficientsP() -> List[float]:
    """Gets all 'p' coefficients

    Returns
    ------- 
    list of float
        A list containing all 'p' coefficients
    """
    return [4.06054, 12.8802, -30.5858, 16.7498]

def getApolloCoefficients() -> List[float]:
    """Coefficients used for the adjustment of the ROLO model using Apollo spectra.

    Returns
    -------
    list of float
        A list containing all Apollo coefficients
    """
    return [1.0301, 1.0970, 0.9325, 0.9466, 1.0225, 1.0157, 1.0470, 1.0084, 1.0100, 1.0148, 0.9843, 1.0134, 0.9329, 0.9849, 0.9994, 0.9957, 1.0059, 0.9618, 0.9561, 0.9796, 0.9568, 0.9873, 1.0575, 1.0108, 0.9743, 1.0386, 1.0338, 1.0577, 1.0650, 1.0815, 0.8945, 0.9689]