"""ESI Extraterrestrial Solar Irradiation

This module contains the functionalities for obtaining the extraterrestrial solar irradiation
of a concrete wavelength, based on Wehrli (1985).

It exports the following classes:
    * WehrliFile - Enum that represents which wehrli data source will be used in the calculation
        of the ESI
    * ESIMethod - Enum that represents which interpolation method will be used in the calculation
        of the ESI
    * GaussianFilterParams - Dataclass that represents the parameters for the gaussian filter
        interpolation
    * ESICalculator - Calculator of Extraterrestrial Solar Irradiance.
"""
import csv
from dataclasses import dataclass
from io import StringIO
from typing import Tuple, Dict, List, Iterable, Optional
from threading import RLock
import math
import pkgutil
import enum
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray


def _linear_interpolation(
    wavelength_nm: float, x_values: List[float], y_values: List[float]
) -> float:
    """
    Wrapper that linearly interpolates

    Parameters
    ----------
    wavelength_nm : float
        x value that is going to be interpolated
    x_values : list of float
        list of x values (keys)
    y_values : list of float
        list of y values

    Returns
    -------
    float
        Interpolated value.
    """
    # This works because, supposedly, python dicts preserve insertion order since 3.7
    return np.interp(wavelength_nm, x_values, y_values)


def _gaussian_filter_non_eq_filter_input(
    center: float, all_x: List[float], all_y: List[float], radius: float, sigma: float
) -> Tuple[List[float], List[float]]:
    """
    Calculates the gaussian values from all_x that are in the range of
    [center-radius, center+radius].

    Parameters
    ----------
    center : float
        Center of the gaussian distribution
    all_x : list of float
        All x values (keys) present on the gaussian distribution
    all_y : list of float
        All y values, in the same order as all_y
    radius : float
        Radius that will filter the accepted values. From all_x, each one that is in the
            interval [center-radius, center+radius]
    sigma : float
        Standard deviation of the Gaussian Filter.

    Returns
    -------
    list of float
        Gaussian values of the elements from all_x in the range.
    list of float
        y values relative to the filtered values, in the same order.
    """
    min_x = center - radius
    max_x = center + radius
    gauss_vals = []
    final_y = []
    sigma = float(sigma)
    for i, x_val in enumerate(all_x):
        if min_x <= x_val <= max_x:
            gauss_param = float(x_val - center)
            val = (1 / (sigma * math.sqrt(2 * math.pi))) * (
                math.exp(-(gauss_param**2) / (2 * sigma**2))
            )
            gauss_vals.append(val)
            final_y.append(all_y[i])
    return gauss_vals, final_y


def _gaussian_filter_non_eq_sum_percentages(
    gauss_vals: List[float], final_y: List[float]
) -> float:
    """
    Adjust the gaussian values for the non equidistant data, and calculates
    the gaussian filtered value.

    Parameters
    ----------
    gauss_vals : list of float
        list of gaussian values used for the calculation of the final value.
    final_y : list of float
        list with the y values relative to the gaussian values.

    Returns
    -------
    float
        Gaussian-filtered value
    """
    gauss_sum = sum(gauss_vals)
    val_sum = 0
    for i, final_y_val in enumerate(final_y):
        if gauss_sum == 0:
            perc = 0
        else:
            perc = gauss_vals[i] / gauss_sum
        val_sum += perc * final_y_val
    return val_sum


def _gaussian_filter_non_equidistant(
    center: float,
    all_x: List[float],
    all_y: List[float],
    radius: float = 1,
    sigma: float = 1,
) -> float:
    """
    Gaussian filter for non equidistant data.

    Parameters
    ----------
    center : float
        Center of the gaussian distribution
    all_x : list of float
        All x values (keys) present on the gaussian distribution
    all_y : list of float
        All y values, in the same order as all_y
    radius : float
        Radius that will filter the accepted values. From all_x, each one that is in the
            interval [center-radius, center+radius]
    sigma : float
        Standard deviation of the Gaussian Filter.

    Returns
    -------
    float
        Gaussian-filtered value
    """
    gauss_vals, final_y = _gaussian_filter_non_eq_filter_input(
        center, all_x, all_y, radius, sigma
    )
    return _gaussian_filter_non_eq_sum_percentages(gauss_vals, final_y)


class WehrliFile(enum.Enum):
    """
    Wehrli data location that will be used in the calculation of the ESI.

    Values
    ------
    ORIGINAL_WEHRLI : Original wehrli data.
    SIMPLE_FILTER_WEHRLI : Wehrli data passed through a gaussian filter and linear interpolation.
        (See utils/wehrli_gauss).
    GAUSSIAN_WEHRLI : Wehrli data passed through a gaussian filter with data for every 0.1 nm.
        Similar effectiveness to using ORIGINAL_WEHRLI with GAUSSIAN_FILTER.
        Not recommended for obtaining "Wm⁻²" data, only "Wm⁻²/nm".
    ASC_WEHRLI : Wehrli data passed through different filters. This is the default one as it's the
        one used in AEMET's RimoApp. Not recommended for obtaining "Wm⁻²" data, only "Wm⁻²/nm".
    """

    ORIGINAL_WEHRLI = "data/wehrli_original.csv"
    SIMPLE_FILTER_WEHRLI = "data/wehrli_filtered.csv"
    GAUSSIAN_WEHRLI = "data/wehrli_gaussian.csv"
    ASC_WEHRLI = "data/wehrli_asc.csv"


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


@dataclass
class GaussianFilterParams:
    """
    Parameters for the gaussian filter interpolation

    Attributes
    ----------
    radius : float
        Radius of the width of the Gaussian filter.
    sigma : float
        Standard deviation for the Gaussian filter.
    """

    radius: float = 1
    sigma: float = 1


class ESICalculator(ABC):
    """
    Abstract Calculator of Extraterrestrial Solar Irradiance.
    """

    @abstractmethod
    def get_esi(
        self, wavelengths_nm: Iterable[float], per_nm: bool = False
    ) -> NDArray[np.float32]:
        """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
        Returns the data in Wm⁻² or Wm⁻²/nm

        Parameters
        ----------
        wavelengths_nm : iterable of float
            Wavelengths (in nanometers) of which the extraterrestrial solar irradiance will be
            obtained
        per_nm : bool
            If True the irradiance will be in Wm⁻²/nm, otherwise it will be in Wm⁻².
            Default is False.

        Returns
        -------
        np.array of float
            The expected extraterrestrial solar irradiance
        """


class ESICalculatorWehrli(ESICalculator):
    """
    Calculator of Extraterrestrial Solar Irradiance.
    Based on Wehrli data and some sort of interpolation.

    Attributes
    ----------
    wehrli_file : WehrliFile
        Wehrli data source that will be used in the calculation of the ESI. It could be the original
        data or some filtered data.
    method : ESIMethod
        Interpolation method that will be used in the calculation of the ESI.
    gfp : GaussianFilterParams
        Parameters of the gaussian filter method, in case that that is the chosen one.
    """

    __slots__ = ["wehrli_file", "method", "gfp"]
    _cache_wehrli_files: Dict[str, Dict[float, Tuple[float, float]]] = {}
    _cache_lock: RLock = RLock()

    def __init__(
        self,
        wehrli_file: WehrliFile = WehrliFile.ASC_WEHRLI,
        method: ESIMethod = ESIMethod.LINEAR_INTERPOLATION,
        gaussian_filter_params: Optional[GaussianFilterParams] = None,
    ):
        """
        Parameters
        ----------
        wehrli_file : WehrliFile
            Wehrli data source that will be used in the calculation of the ESI. It could be the
            original data or some filtered data.
        method : ESIMethod
            Interpolation method that will be used in the calculation of the ESI.
        gfp : GaussianFilterParams
            Parameters of the gaussian filter method, in case that that is the chosen one.
            Default = None.
        """
        self.wehrli_file = wehrli_file
        self.method = method
        if gaussian_filter_params is None:
            self.gfp = GaussianFilterParams()
        else:
            self.gfp = gaussian_filter_params

    @classmethod
    def _store_wehrli_cache(cls, filename: str, data: Dict[float, Tuple[float, float]]):
        with cls._cache_lock:
            cls._cache_wehrli_files[filename] = data

    @classmethod
    def _read_wehrli_cache(
        cls, filename: str
    ) -> Optional[Dict[float, Tuple[float, float]]]:
        with cls._cache_lock:
            return cls._cache_wehrli_files.get(filename)

    @classmethod
    def _clear_cache(cls):
        with cls._cache_lock:
            cls._cache_wehrli_files.clear()

    def _get_wehrli_data(self) -> Dict[float, Tuple[float, float]]:
        """Returns all wehrli data

        Returns
        -------
        A dict that has the wavelengths as keys (float), and as values it has tuples of the
        (Wm⁻²/nm, Wm⁻²) values.
        """
        data = self._read_wehrli_cache(self.wehrli_file.value)
        if data is not None:
            return data
        wehrli_bytes = pkgutil.get_data(__name__, self.wehrli_file.value)
        wehrli_string = wehrli_bytes.decode()
        with StringIO(wehrli_string) as file:
            csvreader = csv.reader(file)
            next(csvreader)  # Discard the header
            data = {}
            for row in csvreader:
                data[float(row[0])] = (float(row[1]), float(row[2]))
        self._store_wehrli_cache(self.wehrli_file.value, data)
        return data

    def get_esi(
        self, wavelengths_nm: Iterable[float], per_nm: bool = False
    ) -> NDArray[np.float32]:
        """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
        Returns the data in Wm⁻² or Wm⁻²/nm

        Parameters
        ----------
        wavelengths_nm : iterable of float
            Wavelengths (in nanometers) of which the extraterrestrial solar irradiance will be
            obtained
        per_nm : bool
            If True the irradiance will be in Wm⁻²/nm, otherwise it will be in Wm⁻².
            Default is False.

        Returns
        -------
        np.array of float
            The expected extraterrestrial solar irradiance
        """
        wavelengths_nm = np.array(wavelengths_nm)
        wehrli_data = self._get_wehrli_data()
        wehrli_x = list(wehrli_data.keys())
        value_index = 1
        if per_nm:
            value_index = 0
        wavelengths_nm = np.where(
            wavelengths_nm < wehrli_x[0], wehrli_x[0], wavelengths_nm
        )
        wavelengths_nm = np.where(
            wavelengths_nm > wehrli_x[-1], wehrli_x[-1], wavelengths_nm
        )
        wehrli_y = list(map(lambda x: x[value_index], wehrli_data.values()))
        if self.method == ESIMethod.LINEAR_INTERPOLATION:
            return _linear_interpolation(wavelengths_nm, wehrli_x, wehrli_y)
        # ESIMethod.GAUSSIAN_FILTER
        gauss_res = _gaussian_filter_non_equidistant(
            wavelengths_nm, wehrli_x, wehrli_y, self.gfp.radius, self.gfp.sigma
        )
        if (
            gauss_res == 0
        ):  # There was no wehrli data near enough from the given wavelength_nm
            return _linear_interpolation(wavelengths_nm, wehrli_x, wehrli_y)
        return gauss_res
