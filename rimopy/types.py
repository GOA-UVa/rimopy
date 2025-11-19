"""types

This module contains the class 'MoonDatas', that conveys the needed data for the calculation of
extraterrestrial lunar irradiance. The data is probably obtained from NASA's SPICE Toolbox

It exports the following classes:
    * MoonDatas - Moon data needed to calculate Moon's irradiance.
    * MissingRCFBehavior - Enum with the options to do when an RCF is missing for a wavelength.
    * EarthPoint - Data of the point on Earth surface of which RIMO will be computed.
"""

from typing import Iterable, Union, List
from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy.typing import NDArray


class MoonDatas:
    """
    Moon data needed to calculate Moon's irradiance, probably obtained from NASA's SPICE Toolbox

    Attributes
    ----------
    distance_sun_moon : array of float
        Distance between the Sun and the Moon (in astronomical units)
    distance_observer_moon : array of float
        Distance between the Observer and the Moon (in kilometers)
    long_sun_radians : array of float
        Selenographic longitude of the Sun (in radians)
    lat_obs : array of float
        Selenographic latitude of the observer (in degrees)
    long_obs : array of float
        Selenographic longitude of the observer (in degrees)
    mpa_degrees: array of float
        Moon phase angle (in degrees)
    absolute_mpa_degrees : array of float
        Absolute Moon phase angle (in degrees)
    """

    def __init__(
        self,
        dsm: Iterable[float],
        dom: Iterable[float],
        lonsun: Iterable[float],
        latobs: Iterable[float],
        lonobs: Iterable[float],
        mpa: Iterable[float],
    ):
        ampa = ((np.abs(mpa) + 180) % 360) - 180
        self._data = np.array([dsm, dom, lonsun, latobs, lonobs, mpa, ampa])

    @property
    def dsm(self) -> NDArray[np.float64]:
        return self._data[0]

    @property
    def dom(self) -> NDArray[np.float64]:
        return self._data[1]

    @property
    def lonsun(self) -> NDArray[np.float64]:
        return self._data[2]

    @property
    def latobs(self) -> NDArray[np.float64]:
        return self._data[3]

    @property
    def lonobs(self) -> NDArray[np.float64]:
        return self._data[4]

    @property
    def mpa(self) -> NDArray[np.float64]:
        return self._data[5]

    @property
    def ampa(self) -> NDArray[np.float64]:
        return self._data[6]

    def get_moondata(self, i) -> NDArray[np.float64]:
        return self._data[:, i]


class MissingRCFBehavior(str, Enum):
    """What to do when an RCF is missing for a wavelength."""

    ERROR = "error"  # raise ValueError
    WARN = "warn"  # issue a warning, return uncorrected values
    IGNORE = "ignore"  # silently return uncorrected values
    NEAREST = "nearest"  # return the RCF of the closest wavelength with an existing one


@dataclass
class EarthPoint:
    """
    Data of the point on Earth surface of which the ELI will be calculated.

    Attributes
    ----------
    lat : float
        Geographic latitude (in degrees) of the location.
    lon : float
        Geographic longitude (in degrees) of the location.
    utc_times : list of str | str
        Time/s at which the ELI will be calculated, in a valid UTC DateTime format.
    altitude : float
        Altitude over the sea level in meters. Default = 0.
    """

    __slots__ = ["lat", "lon", "utc_times", "altitude"]

    def __init__(
        self,
        lat: float,
        lon: float,
        utc_times: Union[List[str], str],
        altitude: float = 0,
    ):
        """
        Parameters
        ----------
        lat : float
            Geographic latitude (in degrees) of the location.
        lon : float
            Geographic longitude (in degrees) of the location.
        utc_times : list of str | str
            Time/s at which the ELI will be calculated, in a valid UTC DateTime format.
        altitude : float
            Altitude over the sea level in meters. Default = 0.
        """
        self.lat = lat
        self.lon = lon
        self.altitude = altitude
        if isinstance(utc_times, list):
            self.utc_times = utc_times
        else:
            self.utc_times = [utc_times]

    def set_utc_times(self, utc_times: Union[List[str], str]):
        """
        Modifies the utc_times attribute

        Parameters
        ----------
        utc_times : list of str | str
            Time/s at which the ELI will be calculated, in a valid UTC DateTime format.
        """
        if isinstance(utc_times, list):
            self.utc_times = utc_times
        else:
            self.utc_times = [utc_times]
