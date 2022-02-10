"""ESI Extraterrestrial Solar Irradiation

This module contains the functionalities for obtaining the extraterrestrial solar irradiation
of a concrete wavelength, based on Wehrli (1985).

It exports the foollowing functions:

    * getESI - returns the expected extraterrestrial solar irradiation of a concrete wavelength in Wm⁻²
	* getESIPerNm - returns the expected extraterrestrial solar irradiation of a concrete wavelength in Wm⁻²/nm
"""

import csv
from io import StringIO
from typing import Tuple, Dict, List
import math
import pkgutil
from scipy.interpolate import interp1d
import enum

def _linearInterpolation(wavelength_nm: float, x: List[float], y: List[float]):
    f = interp1d(x, y, 'linear') # This works because, supposedly, python dicts preserve insertion order since 3.7
    return f(wavelength_nm).item()

def _gaussianFilteredNonEquidistant(center: float, all_x: List[float], all_y: List[float], radius=1, sigma=1):
    min_x = center - radius
    max_x = center + radius
    gauss_vals = []
    final_y = []
    gauss_sum = 0
    for i in range(len(all_x)):
        if all_x[i] >= min_x and all_x[i] <= max_x:
            gauss_param = all_x[i] - center
            val = (1/(sigma*math.sqrt(2*math.pi)))*(math.exp(-(gauss_param**2)/(2*sigma**2)))
            gauss_vals.append(val)
            gauss_sum += val
            final_y.append(all_y[i])
    val_sum = 0
    for i in range(len(final_y)):
        if gauss_sum == 0: perc = 0
        else: perc = gauss_vals[i]/gauss_sum
        val_sum += perc * final_y[i]
    return val_sum

class WehrliFile(enum.Enum):
    """
    Wehrli data location that will be used in the calculation of the ESI.

    Values
    ------
    ORIGINAL_WEHRLI : Original wehrli data.
    SIMPLE_FILTER_WEHRLI : Wehrli data passed through a gaussian filter and linear interpolation. (See utils/wehrli_gauss).
    """
    ORIGINAL_WEHRLI = 'data/wehrli_original.csv'
    SIMPLE_FILTER_WEHRLI = 'data/wehrli_filtered.csv'

class ESIMethod(enum.Enum):
    """
    Interpolation method that will be used in the calculation of the ESI.

    Values
    ------
    LINEAR_INTERPOLATION : The method will be linear interpolation.
    GAUSSIAN_FILTER : The method will be a gaussian filter.
    """
    LINEAR_INTERPOLATION = 1
    GAUSSIAN_FILTER = 2

class GaussianFilterParams():
    """
    Parameters for the gaussian filter interpolation

    Attributes
    ----------
    radius : float
        Radius of the width of the Gaussian filter.
    sigma : float
        Standard deviation for the Gaussian filter.
    """
    def __init__(self, radius: float = 1, sigma: float = 1):
        """
        Parameters
        ----------
        radius : float
            Radius of the width of the Gaussian filter.
        sigma : float
            Standard deviation for the Gaussian filter.
        """
        self.radius = radius
        self.sigma = sigma

class ESICalculator():
    """
    Calculator of Extraterrestrial Solar Irradiance.
    Based on Wehrli data and some sort of interpolation.

    Attributes
    ----------
    wehrli_file : WehrliFile
        Wehrli data source that will be used in the calculation of the ESI. It could be the original data or some filtered data.
    method : ESIMethod
        Interpolation method that will be used in the calculation of the ESI.
    gfp : GaussianFilterParams
        Parameters of the gaussian filter method, in case that that is the chosen one.
    """
    __slots__ = ['wehrli_file', 'method', 'gfp']
    def __init__(self, wehrli_file: WehrliFile=WehrliFile.SIMPLE_FILTER_WEHRLI, method: ESIMethod=ESIMethod.LINEAR_INTERPOLATION, gaussian_filter_params: GaussianFilterParams=None):
        """
        Parameters
        ----------
        wehrli_file : WehrliFile
            Wehrli data source that will be used in the calculation of the ESI. It could be the original data or some filtered data.
        method : ESIMethod
            Interpolation method that will be used in the calculation of the ESI.
        gfp : GaussianFilterParams
            Parameters of the gaussian filter method, in case that that is the chosen one. Default = None.
        """
        self.wehrli_file = wehrli_file
        self.method = method
        if gaussian_filter_params == None:
            self.gfp = GaussianFilterParams()
        else: self.gfp = gaussian_filter_params
    
    def _getWehrliData(self) -> Dict[float, Tuple[float, float]]:
        """Returns all wehrli data

        Returns
        -------
        A dict that has the wavelengths as keys (float), and as values it has tuples of the (Wm⁻²/nm, Wm⁻²/sm) values.
        """
        wehrli_bytes = pkgutil.get_data(__name__, self.wehrli_file.value)
        wehrli_string = wehrli_bytes.decode()
        file = StringIO(wehrli_string)
        csvreader = csv.reader(file)
        next(csvreader) # Discard the header
        data = {}
        for row in csvreader:
            data[float(row[0])] = (float(row[1]), float(row[2]))
        file.close()
        return data

    def getESI(self, wavelength_nm: float) -> float:
        """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
        Returns the data in Wm⁻²

        Parameters
        ----------
        wavelength_nm : float
            Wavelength (in nanometers) of which the extraterrestrial solar irradiance will be obtained

        Returns
        -------
        float
            The expected extraterrestrial solar irradiance in Wm⁻²
        """
        wehrli_data = self._getWehrliData()
        wehrli_x = list(wehrli_data.keys())
        if wavelength_nm in wehrli_x:
            return wehrli_data[wavelength_nm][1]
        if wavelength_nm < wehrli_x[0]:
            return wehrli_data[wehrli_x[0]][1]
        if wavelength_nm > wehrli_x[-1]:
            return wehrli_data[wehrli_x[-1]][1]
        wehrli_y = list(map(lambda x : x[1], wehrli_data.values()))
        if self.method == ESIMethod.LINEAR_INTERPOLATION:
            return _linearInterpolation(wavelength_nm, wehrli_x, wehrli_y)
        else:
            gauss_res =  _gaussianFilteredNonEquidistant(wavelength_nm, wehrli_x, wehrli_y, self.gfp.radius, self.gfp.sigma)
            if gauss_res == 0: # There was no wehrli data near enough from the given wavelength_nm
                return _linearInterpolation(wavelength_nm, wehrli_x, wehrli_y)
            return gauss_res


    def getESIPerNm(self, wavelength_nm: float) -> float:
        """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
        Returns the data in Wm⁻²/nm

        Parameters
        ----------
        wavelength_nm : float
            Wavelength (in nanometers) of which the extraterrestrial solar irradiance will be obtained

        Returns
        -------
        float
            The expected extraterrestrial solar irradiance in Wm⁻²/nm
        """
        wehrli_data = self._getWehrliData()
        wehrli_x = list(wehrli_data.keys())
        if wavelength_nm in wehrli_x:
            return wehrli_data[wavelength_nm][1]
        if wavelength_nm < wehrli_x[0]:
            return wehrli_data[wehrli_x[0]][1]
        if wavelength_nm > wehrli_x[-1]:
            return wehrli_data[wehrli_x[-1]][1]
        wehrli_y = list(map(lambda x : x[0], wehrli_data.values()))
        if self.method == ESIMethod.LINEAR_INTERPOLATION:
            return _linearInterpolation(wavelength_nm, wehrli_x, wehrli_y)
        else:
            gauss_res = _gaussianFilteredNonEquidistant(wavelength_nm, wehrli_x, wehrli_y, self.gfp.radius, self.gfp.sigma)
            if gauss_res == 0: # There was no wehrli data near enough from the given wavelength_nm
                return _linearInterpolation(wavelength_nm, wehrli_x, wehrli_y)
            return gauss_res