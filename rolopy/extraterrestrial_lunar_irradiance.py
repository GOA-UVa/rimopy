"""Extraterrestrial Lunar Irradiance

This module is important
"""

import math
from typing import List
import coefficients as coeffs

class Moon_Data:
    """
    Moon data needed to calculate Moon's irradiance, probably obtained from NASA's SPICE Toolbox

    ...

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
    """
    __slots__ = ['distance_sun_moon', 'distance_observer_moon', 'long_sun_radians', 'lat_obs', 'long_obs']

    def __init__(self, distance_sun_moon: float, distance_observer_moon: float, long_sun_radians: float, lat_obs: float, long_obs: float):
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
        """
        self.distance_sun_moon = distance_sun_moon
        self.distance_observer_moon = distance_observer_moon
        self.long_sun_radians = long_sun_radians
        self.lat_obs = lat_obs
        self.long_obs = long_obs

def summatory_a(k: int, gr: float) -> float:
    count: float = 0.0
    a: List[List[float]] = coeffs.getCoefficientsA()
    for i in range (len(a)):
        count = count + a[i][k] * gr ^ i 
    return count

def summatory_b(k: int, phi: float) -> float:
    count: float = 0.0
    b: List[List[float]] = coeffs.getCoefficientsB()
    for j in range (len(b)):
        count = count + b[j][k] * phi ^ (2*(j + 1) - 1)
    return count

def ln_moon_disk_reflectance(absolute_MPA_degrees: float, wavelength_index: int, moon_data: 'Moon_Data') -> float:
    k = wavelength_index
    gd = absolute_MPA_degrees
    gr = math.radians(gd)
    phi = moon_data.long_sun_radians
    c: List[float] = coeffs.getCoefficientsC()
    d: List[List[float]] = coeffs.getCoefficientsD()
    p: List[float] = coeffs.getCoefficientsP()
    l_theta = moon_data.lat_obs
    l_phi = moon_data.long_obs
    sum_a = summatory_a(k, gr)
    sum_b = summatory_b(k, phi)
    d1 = d[0][k] * math.exp( - gd / p[0])
    d2 = d[1][k] * math.exp( - gd / p[1])
    d3 = d[2][k] * math.cos( (gd - p[2]) / p[3])
    result = sum_a + sum_b + c[0] * l_phi + c[1] * l_theta + c[2] * phi * l_phi + c[3] * phi * l_theta + d1 + d2 + d3
    return result

def getExtraterrestrialSolarIrradiance(wavelength_nm: float) -> float:
    pass

def getIrradianceForWavelength(wavelength_nm: float, absolute_MPA_degrees: float, moon_data: 'Moon_Data') -> float:
    index = coeffs.getWavelengths().index(wavelength_nm)
    lnAk = ln_moon_disk_reflectance(absolute_MPA_degrees, index, moon_data)
    ak = math.exp(lnAk)
    solid_angle_moon: float = 6.4177 * 10 ^ -5
    distance_earth_moon_km: int = 384400
    omega = solid_angle_moon
    esk = getExtraterrestrialSolarIrradiance(wavelength_nm)
    dsm = moon_data.distance_sun_moon
    dom = moon_data.distance_observer_moon
    em = ((ak * omega * esk) / math.pi) * ((1 / dsm) ^ 2) * (distance_earth_moon_km / dom) ^ 2
    return em