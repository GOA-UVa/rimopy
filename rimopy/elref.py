"""ELRef Extraterrestrial Lunar Reflectance

This module allows for the calculation of the extraterrestrial lunar reflectance.
"""

from typing import Iterable

import numpy as np

from . import coefficients as coeffs
from .types import MoonData
from . import coefficients as coeffs
from . import correction_factor as corr_f


def _summatory_a(wavelength_nm: float, gr_value: float) -> float:
    """The first summatory of Eq. 2 in Roman et al., 2020

    Parameters
    ----------
    wavelength_nm : float
        Wavelength in nanometers from which the moon's disk reflectance is being calculated
    gr_value : float
        Absolute value of MPA in radians

    Returns
    -------
    float
        Result of the computation of the first summatory
    """
    count: float = 0.0
    a_coeffs: Iterable[float] = coeffs.get_coefficients_a(wavelength_nm)
    for i, a_value in enumerate(a_coeffs):
        count = count + a_value * gr_value ** i
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
    b_coeffs: Iterable[float] = coeffs.get_coefficients_b(wavelength_nm)
    for j, b_value in enumerate(b_coeffs):
        count = count + b_value * phi ** (2*(j + 1) - 1)
    return count


def _get_correction_factor(wavelength_nm: float, mpa: float) -> float:
    """Calculation of RIMO correction factor (RCF) following Eq 9 in Roman et al., 2020

    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated
    mpa : float
        Absolute Moon phase angle (in radians)

    Returns
    -------
    float
        The calculated RCF
    """
    params = corr_f.get_correction_params(wavelength_nm)
    rcf = params.a_coeff + params.b_coeff *mpa + params.c_coeff * mpa ** 2
    return rcf


def _ln_moon_disk_reflectance(absolute_mpa_degrees: float, wavelength_nm: float,
                              moon_data: MoonData) -> float:
    """The calculation of the ln of the reflectance of the Moon's disk, following Eq.2 in
    Roman et al., 2020

    If the wavelength has no associated ROLO coefficients, it uses some linearly interpolated
    ones.

    Parameters
    ----------
    absolute_mpa_degrees : float
        Absolute Moon phase angle (in degrees)
    wavelength_nm : float
        Wavelength in nanometers from which one wants to obtain the MDR.
    moon_data : MoonData
        Moon data needed to calculate Moon's irradiance

    Returns
    -------
    float
        The ln of the reflectance of the Moon's disk for the inputed data
    """
    gd_value = absolute_mpa_degrees
    gr_value = np.radians(gd_value)
    phi = moon_data.long_sun_radians
    c_coeffs: Iterable[float] = coeffs.get_coefficients_c()
    d_coeffs: Iterable[float] = coeffs.get_coefficients_d(wavelength_nm)
    p_coeffs: Iterable[float] = coeffs.get_coefficients_p()
    l_theta = moon_data.lat_obs
    l_phi = moon_data.long_obs
    sum_a = _summatory_a(wavelength_nm, gr_value)
    sum_b = _summatory_b(wavelength_nm, phi)
    d1_value = d_coeffs[0] * np.exp(- gd_value / p_coeffs[0])
    d2_value = d_coeffs[1] * np.exp(- gd_value / p_coeffs[1])
    d3_value = d_coeffs[2] * np.cos((gd_value - p_coeffs[2]) / p_coeffs[3])
    result = sum_a + sum_b + c_coeffs[0] * l_phi + c_coeffs[1] * l_theta + c_coeffs[2] * phi * \
        l_phi + c_coeffs[3] * phi * l_theta + d1_value + d2_value + d3_value
    return result


def _interpolated_moon_disk_reflectance(absolute_mpa_degrees: float, wavelength_nm: float,
                                        moon_data: 'MoonData', adjust_apollo: bool) -> float:
    """The calculation of the reflectance of the Moon's disk, following Eq.2 in Roman et al., 2020

    If the wavelength is not present in the ROLO coefficients, it calculates the linear
    interpolation between the previous and the next one, or the extrapolation with the two
    nearest ones in case that it's on an extreme.

    Parameters
    ----------
    absolute_mpa_degrees : float
        Absolute Moon phase angle (in degrees)
    wavelength_nm : float
        Wavelength in nanometers from which one wants to obtain the MDR.
    moon_data : 'MoonData'
        Moon data needed to calculate Moon's irradiance
    adjust_apollo : bool
        If True, the calculated reflectance will be adjusted to the Apollo spectra.

    Returns
    -------
    float
        The ln of the reflectance of the Moon's disk for the inputed data
    """
    wvlens = coeffs.get_wavelengths()
    if adjust_apollo:
        apollo_coeffs = coeffs.get_apollo_coefficients()
    else:
        apollo_coeffs = [1 for i in range(len(wvlens))]
    if wavelength_nm in wvlens:
        apollo_i = wvlens.index(wavelength_nm)
        return np.exp(_ln_moon_disk_reflectance(absolute_mpa_degrees, wavelength_nm,
                                                  moon_data))*apollo_coeffs[apollo_i]
    if wavelength_nm < wvlens[0]:
        # The extrapolation done is "nearest"
        return _interpolated_moon_disk_reflectance(absolute_mpa_degrees, wvlens[0],
                                                   moon_data, adjust_apollo)
    if wavelength_nm > wvlens[-1]:
        # The extrapolation done is "nearest"
        return _interpolated_moon_disk_reflectance(absolute_mpa_degrees, wvlens[-1],
                                                   moon_data, adjust_apollo)
    near_left = -np.inf
    near_right = np.inf
    for wvlen in wvlens:
        if near_left < wvlen < wavelength_nm:
            near_left = wvlen
        elif  wavelength_nm < wvlen < near_right:
            near_right = wvlen
    x_values = [near_left, near_right]
    left_index = wvlens.index(x_values[0])
    right_index = wvlens.index(x_values[1])
    y_values = []
    y_values.append(np.exp(_ln_moon_disk_reflectance(absolute_mpa_degrees, x_values[0],
                                                       moon_data)) * apollo_coeffs[left_index])
    y_values.append(np.exp(_ln_moon_disk_reflectance(absolute_mpa_degrees, x_values[1],
                                                       moon_data)) * apollo_coeffs[right_index])
    return np.interp(wavelength_nm, x_values, y_values)


def get_reflectance_interpolating_coefficients(absolute_mpa_degrees: float, wavelength_nm: float, moon_data: MoonData, apply_correction: bool = True):
    a_l = np.exp(_ln_moon_disk_reflectance(absolute_mpa_degrees, wavelength_nm, moon_data))
    if apply_correction:
        mr_correction_factor = _get_correction_factor(
            wavelength_nm,
            np.radians(moon_data.absolute_mpa_degrees)
        )
        a_l = a_l * mr_correction_factor
    return a_l


def get_interpolated_reflectance(absolute_mpa_degrees: float, wavelength_nm: float, moon_data: MoonData, apply_correction: bool = True, adjust_apollo: bool = True):
    a_l = _interpolated_moon_disk_reflectance(absolute_mpa_degrees, wavelength_nm, moon_data, adjust_apollo)
    if apply_correction:
        mr_correction_factor = _get_correction_factor(
            wavelength_nm,
            np.radians(moon_data.absolute_mpa_degrees)
        )
        a_l = a_l * mr_correction_factor
    return a_l
