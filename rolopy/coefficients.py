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
"""
from typing import List

def getWavelengths() -> List[float]:
    """Gets all wavelengths present in the model, in nanometers

    Returns
    -------
    list of float
        a list of floats that are the wavelengths in nanometers, in order
    """
    return [350.0, 355.1, 405.0, 412.3, 414.4, 441.6, 465.8, 475.0, 486.9]

def getAllCoefficientsA() -> List[List[float]]:
    """Gets all 'a' coefficients

    Returns
    ------- 
    list of list of float
        a list containing multiple list of floats. Each sublist is the list of 'a' coefficients for a wavelength
    """
    return [ [-2.67511, -1.78539, 0.50612, -0.25578], [-2.71924, -1.74298, 0.44523, -0.23315], [-2.35754, -1.72134, 0.40337, -0.21105],
        [-2.34185, -1.74337, 0.42156, -0.21512], [-2.43367, -1.72184, 0.43600, -0.22675], [-2.31964, -1.72114, 0.37286, -0.19304],
        [-2.35085, -1.66538, 0.41802, -0.22541], [-2.28999, -1.63180, 0.36193, -0.20381], [-2.23351, -1.68573, 0.37632, -0.19877 ] ]

def getAllCoefficientsB() -> List[List[float]]:
    """Gets all 'b' coefficients

    Returns
    ------- 
    list of list of float
        a list containing multiple list of floats. Each sublist is the list of 'b' coefficients for a wavelength
    """
    return [ [0.03744, 0.00981, -0.00322], [0.03492, 0.01142, -0.00383], [0.03505, 0.01043, -0.00341], [0.03141, 0.01364, -0.00472],
        [0.03474, 0.01188, -0.00422 ], [0.03736, 0.01545, -0.00559], [0.04274, 0.01127, -0.00439], [0.04007, 0.01216, -0.00437],
        [0.03881, 0.01566, -0.00555] ]

def getAllCoefficientsD() -> List[List[float]]:
    """Gets all 'd' coefficients

    Returns
    ------- 
    list of list of float
        a list containing multiple list of floats. Each sublist is the list of 'd' coefficients for a wavelength
    """
    return [ [0.34185, 0.01441, -0.01602], [0.33875, 0.01612, -0.00996], [0.35235, -0.03818, -0.00006], [0.36591, -0.05902, 0.00080],
        [0.35558, -0.03247, -0.00503], [0.37935, -0.09562, 0.00970], [0.33450, -0.02546, -0.00484], [0.33024, -0.03131, 0.00222],
        [0.36590, -0.08945, 0.00678] ]

def getCoefficientsA(wavelength_index: int) -> List[float]:
    """Gets all 'a' coefficients for a concrete wavelength

    Parameters
    ----------
    wavelength_index : int
        Index of the wavelength from which one wants to obtain the coefficients. The index is the position of the wavelength on getWavelengths() returned list

    Returns
    ------- 
    list of float
        a list containing the 'a' coefficients for the wavelength
    """
    return getAllCoefficientsA()[wavelength_index]

def getCoefficientsB(wavelength_index: int) -> List[float]:
    """Gets all 'b' coefficients for a concrete wavelength

    Parameters
    ----------
    wavelength_index : int
        Index of the wavelength from which one wants to obtain the coefficients. The index is the position of the wavelength on getWavelengths() returned list

    Returns
    ------- 
    list of float
        a list containing the 'b' coefficients for the wavelength
    """
    return getAllCoefficientsB()[wavelength_index]

def getCoefficientsD(wavelength_index: int) -> List[float]:
    """Gets all 'd' coefficients for a concrete wavelength

    Parameters
    ----------
    wavelength_index : int
        Index of the wavelength from which one wants to obtain the coefficients. The index is the position of the wavelength on getWavelengths() returned list

    Returns
    ------- 
    list of float
        a list containing the 'd' coefficients for the wavelength
    """
    return getAllCoefficientsD()[wavelength_index]

def getCoefficientsC() -> List[float]:
    """Gets all 'c' coefficients

    Returns
    ------- 
    list of float
        a list containing all 'c' coefficients
    """
    return [0.00034115, -0.0013425, 0.00095906, 0.00066229]

def getCoefficientsP() -> List[float]:
    """Gets all 'p' coefficients

    Returns
    ------- 
    list of float
        a list containing all 'p' coefficients
    """
    return [4.06054, 12.8802, -30.5858, 16.7498]