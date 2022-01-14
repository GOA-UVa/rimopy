"""Extraterrestrial Lunar Irradiance

This module is important
"""

import math
from typing import List
from . import coefficients as coeffs
from . import correction_factor as corr_f
from . import esi

class Moon_Data:
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

def summatory_a(wavelength_nm: float, gr: float) -> float:
    """The first summatory of Eq. 2 in Roman et al., 2020

    Parameters
    ----------
    wavelength_nm : float
        Wavelength in nanometers from which the moon's disk reflectance is being calculated
    gr : float
        Absolute value of MPA in radians

    Returns
    -------
    float
        Result of the computation of the first summatory
    """
    count: float = 0.0
    a: List[List[float]] = coeffs.getCoefficientsA(wavelength_nm)
    for i in range (len(a)):
        count = count + a[i] * gr ** i 
    return count

def summatory_b(wavelength_nm: float, phi: float) -> float:
    """The second summatory of Eq. 2 in Roman et al., 2020, without the erratum

    Parameters
    ----------
    wavelength_nm : float
        Wavelength from which the moon's disk reflectance is being calculated
    phi : float
        Selenographic longitude of the Sun (in radians)

    Returns
    -------
    float
        Result of the computation of the second summatory
    """
    count: float = 0.0
    b: List[List[float]] = coeffs.getCoefficientsB(wavelength_nm)
    for j in range (len(b)):
        count = count + b[j] * phi ** (2*(j + 1) - 1)
    return count

def ln_moon_disk_reflectance(absolute_MPA_degrees: float, wavelength_nm: float, moon_data: 'Moon_Data') -> float:
    """The calculation of the ln of the reflectance of the Moon's disk, following Eq.2 in Roman et al., 2020

    Parameters
    ----------
    absolute_MPA_degrees : float
        Absolute Moon phase angle (in degrees)
    wavelength_nm : float
        Wavelength in nanometers from which one wants to obtain the MDR.
    moon_data : 'Moon_Data'
        Moon data needed to calculate Moon's irradiance

    Returns
    -------
    float
        The ln of the reflectance of the Moon's disk for the inputed data
    """
    gd = absolute_MPA_degrees
    gr = math.radians(gd)
    phi = moon_data.long_sun_radians
    c: List[float] = coeffs.getCoefficientsC()
    d: List[float] = coeffs.getCoefficientsD(wavelength_nm)
    p: List[float] = coeffs.getCoefficientsP()
    l_theta = moon_data.lat_obs
    l_phi = moon_data.long_obs
    sum_a = summatory_a(wavelength_nm, gr)
    sum_b = summatory_b(wavelength_nm, phi)
    d1 = d[0] * math.exp( - gd / p[0])
    d2 = d[1] * math.exp( - gd / p[1])
    d3 = d[2] * math.cos( (gd - p[2]) / p[3])
    result = sum_a + sum_b + c[0] * l_phi + c[1] * l_theta + c[2] * phi * l_phi + c[3] * phi * l_theta + d1 + d2 + d3
    return result

def getCorrectionFactor(wavelength_nm: float, mpa: float) -> float:
    """Calculation of RIMO correction factor (RCF) following Eq 9 in Roman et al., 2020

    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial lunar irradiance will be calculated
    mpa : float
        Absolute Moon phase angle (in degrees)

    Returns
    -------
    float
        The calculated RCF
    """
    params = corr_f.getCorrectionParams(wavelength_nm)
    rcf = params.a + params.b *mpa + params.c * mpa ** 2
    return rcf

def getExtraterrestrialSolarIrradiance(wavelength_nm: float) -> float:
    """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
    
    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial solar irradiance will be obtained
    
    Returns
    -------
    float
        The expected extraterrestrial solar irradiance in W/sm
    """
    return esi.getESI(wavelength_nm)

def getIrradianceForWavelength(wavelength_nm: float, absolute_MPA_degrees: float, moon_data: 'Moon_Data') -> float:
    """Calculation of Extraterrestrial Lunar Irradiance following Eq 3 in Roman et al., 2020

    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial lunar irradiance will be calculated
    absolute_MPA_degrees : float
        Absolute Moon phase angle (in degrees)
    moon_data : 'Moon_Data'
        Moon data needed to calculate Moon's irradiance

    Returns
    -------
    float
        The extraterrestrial lunar irradiance calculated
    """
    ln_moon_reflectance = ln_moon_disk_reflectance(absolute_MPA_degrees, wavelength_nm, moon_data)
    mr_correction_factor = getCorrectionFactor(wavelength_nm, absolute_MPA_degrees)
    solid_angle_moon: float = 6.4177 * 10 ** -5

    a_l = math.exp(ln_moon_reflectance) * mr_correction_factor
    omega = solid_angle_moon
    esk = getExtraterrestrialSolarIrradiance(wavelength_nm)
    dsm = moon_data.distance_sun_moon
    dom = moon_data.distance_observer_moon
    distance_earth_moon_km: int = 384400

    em = ((a_l * omega * esk) / math.pi) * ((1 / dsm) ** 2) * (distance_earth_moon_km / dom) ** 2
    return em