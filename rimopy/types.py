"""types

This module contains the class 'MoonDatas', that conveys the needed data for the calculation of
extraterrestrial lunar irradiance. The data is probably obtained from NASA's SPICE Toolbox

It exports the following classes:
    * MoonDatas - Moon data needed to calculate Moon's irradiance.
"""

from typing import Iterable

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

    def __init__(self, dsm: Iterable[float], dom: Iterable[float], lonsun: Iterable[float], latobs: Iterable[float], lonobs: Iterable[float], mpa: Iterable[float]):
        ampa = ((np.abs(mpa)+180) % 360) - 180
        self._data = np.array([dsm, dom, lonsun, latobs, lonobs, mpa, ampa])

    @property
    def dsm(self) -> NDArray[np.float32]:
        return self._data[0]

    @property
    def dom(self) -> NDArray[np.float32]:
        return self._data[1]

    @property
    def lonsun(self) -> NDArray[np.float32]:
        return self._data[2]

    @property
    def latobs(self) -> NDArray[np.float32]:
        return self._data[3]

    @property
    def lonobs(self) -> NDArray[np.float32]:
        return self._data[4]

    @property
    def mpa(self) -> NDArray[np.float32]:
        return self._data[5]

    @property
    def ampa(self) -> NDArray[np.float32]:
        return self._data[6]

    def get_moondata(self, i) -> NDArray[np.float32]:
        return self._data[:, i]
