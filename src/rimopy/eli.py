"""ELI Extraterrestrial Lunar Irradiance

This module is the main module, as it allowes the user to calculate the Extraterrestrial Lunar Irradiance
at a concrete wavelength, at an absolute Moon phase angle, and giving selenographic parameters.

It exports the following functions:

    * getELIBypass - returns the expected extraterrestrial lunar irradiation of a wavelength for any
        observer/solar selenographic coordinates, in Wm⁻².
    * getELI - returns the expected extraterrestrial lunar irradiation of a wavelength in any
        geographic coordinates, in Wm⁻².
    * getELIBypassPerNm - returns the expected extraterrestrial lunar irradiation of a wavelength for any
        observer/solar selenographic coordinates, in Wm⁻²/nm.
    * getELIPerNm - returns the expected extraterrestrial lunar irradiation of a wavelength in any
        geographic coordinates, in Wm⁻²/nm.
"""

import math
from typing import List, Union

from . import spice_iface
from . import coefficients as coeffs
from . import correction_factor as corr_f
from . import esi
from .MoonData import MoonData

def _summatory_a(wavelength_nm: float, gr: float) -> float:
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
    a: List[float] = coeffs.getCoefficientsA(wavelength_nm)
    for i in range (len(a)):
        count = count + a[i] * gr ** i 
    return count

def _summatory_b(wavelength_nm: float, phi: float) -> float:
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
    b: List[float] = coeffs.getCoefficientsB(wavelength_nm)
    for j in range (len(b)):
        count = count + b[j] * phi ** (2*(j + 1) - 1)
    return count

def _ln_moon_disk_reflectance(absolute_MPA_degrees: float, wavelength_nm: float, moon_data: 'MoonData') -> float:
    """The calculation of the ln of the reflectance of the Moon's disk, following Eq.2 in Roman et al., 2020

    Parameters
    ----------
    absolute_MPA_degrees : float
        Absolute Moon phase angle (in degrees)
    wavelength_nm : float
        Wavelength in nanometers from which one wants to obtain the MDR.
    moon_data : 'MoonData'
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
    sum_a = _summatory_a(wavelength_nm, gr)
    sum_b = _summatory_b(wavelength_nm, phi)
    d1 = d[0] * math.exp( - gd / p[0])
    d2 = d[1] * math.exp( - gd / p[1])
    d3 = d[2] * math.cos( (gd - p[2]) / p[3])
    result = sum_a + sum_b + c[0] * l_phi + c[1] * l_theta + c[2] * phi * l_phi + c[3] * phi * l_theta + d1 + d2 + d3
    return result

def _getCorrectionFactor(wavelength_nm: float, mpa: float) -> float:
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

def _getESI(wavelength_nm: float) -> float:
    """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
    Returns the data in Wm⁻²
    
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

def _getESIPerNm(wavelength_nm: float) -> float:
    """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
    Returns the data in Wm⁻²/nm
    
    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial solar irradiance will be obtained
    
    Returns
    -------
    float
        The expected extraterrestrial solar irradiance in W/sm
    """
    return esi.getESIPerNm(wavelength_nm)

def _calculateELI(wavelength_nm: float, moon_data: 'MoonData', applyCorrection: bool, perNm: bool = False) -> float:
    """Calculation of Extraterrestrial Lunar Irradiance following Eq 3 in Roman et al., 2020

    Simulates a lunar observation for a wavelength for any observer/solar selenographic
    latitude and longitude.

    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial lunar irradiance will be calculated
    moon_data : 'MoonData'
        Moon data needed to calculate Moon's irradiance
    applyCorrection : bool
        If True the result will have been multiplied by the RCF (Rimo Correction Factor). Otherwise it won't.
    perNm : bool
        True if the user wants the ELI in Wm⁻²/nm, otherwise it will be in Wm⁻². Default is False.

    Returns
    -------
    float
        The extraterrestrial lunar irradiance calculated
    """
    ln_moon_reflectance = _ln_moon_disk_reflectance(moon_data.absolute_MPA_degrees, wavelength_nm, moon_data)
    solid_angle_moon: float = 6.4177e-05

    a_l = math.exp(ln_moon_reflectance)
    if applyCorrection:
        mr_correction_factor = _getCorrectionFactor(wavelength_nm, moon_data.absolute_MPA_degrees)
        a_l = a_l * mr_correction_factor
    omega = solid_angle_moon
    if perNm:
        esk = _getESIPerNm(wavelength_nm)
    else:
        esk = _getESI(wavelength_nm)
    dsm = moon_data.distance_sun_moon
    dom = moon_data.distance_observer_moon
    distance_earth_moon_km: int = 384400

    em = ((a_l * omega * esk) / math.pi) * ((1 / dsm) ** 2) * (distance_earth_moon_km / dom) ** 2
    return em

def getELIBypass(wavelength_nm: Union[float, List[float]], moon_data: 'MoonData', applyCorrection: bool = True) -> Union[float, List[float]]:
    """Calculation of Extraterrestrial Lunar Irradiance following Eq 3 in Roman et al., 2020

    Allow users to simulate lunar observation for any observer/solar selenographic
    latitude and longitude (thus bypassing the need for their computation from the
    position/time of the observer).

    Returns the data in Wm⁻²

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be calculated.
    moon_data : 'MoonData'
        Moon data needed to calculate Moon's irradiance
    applyCorrection : bool
        If True the result will have been multiplied by the RCF (Rimo Correction Factor). Otherwise it won't.

    Returns
    -------
    float | list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter "wavelength_nm" was a list.
    """
    if isinstance(wavelength_nm, list):
        elis = []
        for w in wavelength_nm:
            elis.append(_calculateELI(w, moon_data, applyCorrection))
        return elis
    return _calculateELI(wavelength_nm, moon_data, applyCorrection)

def getELI(wavelength_nm: Union[float, List[float]], lat: float, long: float, utc_time: str, kernels_path: str, altitude: float = 0, applyCorrection: bool = True) -> Union[float, List[float]]:
    """Calculation of Extraterrestrial Lunar Irradiance from geographic coordinates

    Allow users to simulate lunar observations for any observer position around the Earth
    and at any time.

    Returns the data in Wm⁻²

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be calculated.
    lat : float
        Geographic latitude (in degrees) of the location.
    long : float
        Geographic longitude (in degrees) of the location.
    utc_time : str
        Time at which the ELI will be calculated, in a valid UTC DateTime format.
    kernels_path : str
        Folder where the needed SPICE kernels are stored.
    altitude : float
        Altitude over the sea level in meters. Default = 0.
    applyCorrection : bool
        If True the result will have been multiplied by the RCF (Rimo Correction Factor). Otherwise it won't.

    Returns
    -------
    float | list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter "wavelength_nm" was a list.
    """
    moon_data = spice_iface.getMoonData(lat, long, altitude, utc_time, kernels_path)
    return getELIBypass(wavelength_nm, moon_data, applyCorrection)

def getELIBypassPerNm(wavelength_nm: Union[float, List[float]], moon_data: 'MoonData', applyCorrection: bool = True) -> Union[float, List[float]]:
    """Calculation of Extraterrestrial Lunar Irradiance following Eq 3 in Roman et al., 2020

    Allow users to simulate lunar observation for any observer/solar selenographic
    latitude and longitude (thus bypassing the need for their computation from the
    position/time of the observer).

    Returns the data in Wm⁻²/nm

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be calculated.
    moon_data : 'MoonData'
        Moon data needed to calculate Moon's irradiance
    applyCorrection : bool
        If True the result will have been multiplied by the RCF (Rimo Correction Factor). Otherwise it won't.

    Returns
    -------
    float | list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter "wavelength_nm" was a list.
    """
    if isinstance(wavelength_nm, list):
        elis = []
        for w in wavelength_nm:
            elis.append(_calculateELI(w, moon_data, applyCorrection, True))
        return elis
    return _calculateELI(wavelength_nm, moon_data, applyCorrection, True)

def getELIPerNm(wavelength_nm: Union[float, List[float]], lat: float, long: float, utc_time: str, kernels_path: str, altitude: float = 0, applyCorrection: bool = True) -> Union[float, List[float]]:
    """Calculation of Extraterrestrial Lunar Irradiance from geographic coordinates

    Allow users to simulate lunar observations for any observer position around the Earth
    and at any time.

    Returns the data in Wm⁻²/nm

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be calculated.
    lat : float
        Geographic latitude (in degrees) of the location.
    long : float
        Geographic longitude (in degrees) of the location.
    utc_time : str
        Time at which the ELI will be calculated, in a valid UTC DateTime format.
    kernels_path : str
        Folder where the needed SPICE kernels are stored.
    altitude : float
        Altitude over the sea level in meters. Default = 0.
    applyCorrection : bool
        If True the result will have been multiplied by the RCF (Rimo Correction Factor). Otherwise it won't.

    Returns
    -------
    float | list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter "wavelength_nm" was a list.
    """
    moon_data = spice_iface.getMoonData(lat, long, altitude, utc_time, kernels_path)
    return getELIBypassPerNm(wavelength_nm, moon_data, applyCorrection)

    