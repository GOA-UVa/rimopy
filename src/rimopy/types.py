"""MoonData

This module contains the class 'MoonData', that conveys the needed data for the calculation of
extraterrestrial lunar irradiance. The data is probably obtained from NASA's SPICE Toolbox

It exports the following classes:
    * MoonData - Moon data needed to calculate Moon's irradiance.
"""

from dataclasses import dataclass

@dataclass
class MoonData:
    """
    Moon data needed to calculate Moon's irradiance, probably obtained from NASA's SPICE Toolbox

    Attributes
    ----------
    distance_sun_moon : float
        Distance between the Sun and the Moon (in astronomical units)
    distance_observer_moon : float
        Distance between the Observer and the Moon (in kilometers)
    long_sun_radians : float
        Selenographic longitude of the Sun (in radians)
    lat_obs : float
        Selenographic latitude of the observer (in degrees)
    long_obs : float
        Selenographic longitude of the observer (in degrees)
    absolute_mpa_degrees : float
        Absolute Moon phase angle (in degrees)
    """
    __slots__ = ['distance_sun_moon', 'distance_observer_moon', 'long_sun_radians', 'lat_obs',
        'long_obs', 'absolute_mpa_degrees']

    def __init__(self, distance_sun_moon: float, distance_observer_moon: float,
            long_sun_radians: float, lat_obs: float, long_obs: float, absolute_mpa_degrees: float):
        """
        Parameters
        ----------
        distance_sun_moon : float
            Distance between the Sun and the Moon (in astronomical units)
        distance_observer_moon : float
            Distance between the Observer and the Moon (in kilometers)
        long_sun_radians : float
            Selenographic longitude of the Sun (in radians)
        lat_obs : float
            Selenographic latitude of the observer (in degrees)
        long_obs : float
            Selenographic longitude of the observer (in degrees)
        absolute_mpa_degrees : float
            Absolute Moon phase angle (in degrees)
        """
        self.distance_sun_moon = distance_sun_moon
        self.distance_observer_moon = distance_observer_moon
        self.long_sun_radians = long_sun_radians
        self.lat_obs = lat_obs
        self.long_obs = long_obs
        self.absolute_mpa_degrees = absolute_mpa_degrees
