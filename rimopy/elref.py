"""ELRef Extraterrestrial Lunar Reflectance

This module allows for the calculation of the extraterrestrial lunar reflectance.
"""

from typing import Iterable

import numpy as np
from numpy.typing import NDArray

from . import coefficients as coeffs
from .types import MoonDatas
from . import coefficients as coeffs
from . import correction_factor as corr_f


def _summatory_a(wavelength_nm: float, gr: NDArray[np.float32]) -> NDArray[np.float32]:
    """The first summatory of Eq. 2 in Roman et al., 2020

    Parameters
    ----------
    wavelength_nm : float
        Wavelength in nanometers from which the moon's disk reflectance is being calculated
    gr : array of float
        Absolute value of MPA in radians

    Returns
    -------
    array of float
        Result of the computation of the first summatory
    """
    ac: Iterable[float] = coeffs.get_coefficients_a(wavelength_nm)
    sa = ac[0] + ac[1] * gr + ac[2] * gr**2 + ac[3] * gr**3
    return sa

def _summatory_b(wavelength_nm: float, phi: NDArray[np.float32]) -> NDArray[np.float32]:
    """The second summatory of Eq. 2 in Roman et al., 2020, without the erratum

    Parameters
    ----------
    wavelength_nm : float
        Wavelength from which the moon's disk reflectance is being calculated
    phi : array of float
        Selenographic longitude of the Sun (in radians)

    Returns
    -------
    array of float
        Result of the computation of the second summatory
    """
    bc: Iterable[float] = coeffs.get_coefficients_b(wavelength_nm)
    sb = bc[0] * phi + bc[1] * phi**3 + bc[2] * phi**5
    return sb


def _get_correction_factor(wavelength_nm: float, mpa: NDArray[np.float32]) -> NDArray[np.float32]:
    """Calculation of RIMO correction factor (RCF) following Eq 9 in Roman et al., 2020

    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated
    mpa : array of float
        Absolute Moon phase angle (in radians)

    Returns
    -------
    array of float
        The calculated RCF
    """
    params = corr_f.get_correction_params(wavelength_nm)
    rcf = params.a_coeff + params.b_coeff *mpa + params.c_coeff * mpa ** 2
    return rcf


def _ln_moon_disk_reflectance(absolute_mpa_degrees: NDArray[np.float32], wavelength_nm: float,
                              mds: MoonDatas) -> NDArray[np.float32]:
    """The calculation of the ln of the reflectance of the Moon's disk, following Eq.2 in
    Roman et al., 2020

    If the wavelength has no associated ROLO coefficients, it uses some linearly interpolated
    ones.

    Parameters
    ----------
    absolute_mpa_degrees : array of float
        Absolute Moon phase angle (in degrees)
    wavelength_nm : float
        Wavelength in nanometers from which one wants to obtain the MDR.
    mds : MoonDatas
        Moon data needed to calculate Moon's irradiance

    Returns
    -------
    array of float
        The ln of the reflectance of the Moon's disk for the inputed data
    """
    gd_value = absolute_mpa_degrees
    gr_value = np.radians(gd_value)
    phi = mds.lonsun
    c_coeffs = coeffs.get_coefficients_c()
    d_coeffs = coeffs.get_coefficients_d(wavelength_nm)
    p_coeffs = coeffs.get_coefficients_p()
    l_theta = mds.latobs
    l_phi = mds.lonobs
    sum_a = _summatory_a(wavelength_nm, gr_value)
    sum_b = _summatory_b(wavelength_nm, phi)
    d1_value = d_coeffs[0] * np.exp(- gd_value / p_coeffs[0])
    d2_value = d_coeffs[1] * np.exp(- gd_value / p_coeffs[1])
    d3_value = d_coeffs[2] * np.cos((gd_value - p_coeffs[2]) / p_coeffs[3])
    result = sum_a + sum_b + c_coeffs[0] * l_phi + c_coeffs[1] * l_theta + c_coeffs[2] * phi * \
        l_phi + c_coeffs[3] * phi * l_theta + d1_value + d2_value + d3_value
    return result


def _interpolated_moon_disk_reflectance(absolute_mpa_degrees: NDArray[np.float32], wavelength_nm: float,
                                        mds: MoonDatas, adjust_apollo: bool) -> NDArray[np.float32]:
    """The calculation of the reflectance of the Moon's disk, following Eq.2 in Roman et al., 2020

    If the wavelength is not present in the ROLO coefficients, it calculates the linear
    interpolation between the previous and the next one, or the extrapolation with the two
    nearest ones in case that it's on an extreme.

    Parameters
    ----------
    absolute_mpa_degrees : array of float
        Absolute Moon phase angle (in degrees)
    wavelength_nm : float
        Wavelength in nanometers from which one wants to obtain the MDR.
    mds : MoonDatas
        Moon data needed to calculate Moon's irradiance
    adjust_apollo : bool
        If True, the calculated reflectance will be adjusted to the Apollo spectra.

    Returns
    -------
    array of float
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
                                                  mds))*apollo_coeffs[apollo_i]
    if wavelength_nm < wvlens[0]:
        # The extrapolation done is "nearest"
        return _interpolated_moon_disk_reflectance(absolute_mpa_degrees, wvlens[0],
                                                   mds, adjust_apollo)
    if wavelength_nm > wvlens[-1]:
        # The extrapolation done is "nearest"
        return _interpolated_moon_disk_reflectance(absolute_mpa_degrees, wvlens[-1],
                                                   mds, adjust_apollo)
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
                                                       mds)) * apollo_coeffs[left_index])
    y_values.append(np.exp(_ln_moon_disk_reflectance(absolute_mpa_degrees, x_values[1],
                                                       mds)) * apollo_coeffs[right_index])
    # TODO: Can this be faster?:
    y_values = np.array(y_values).T
    return np.array([np.interp(wavelength_nm, x_values, yval) for yval in y_values])


def get_reflectance_interpolating_coefficients(absolute_mpa_degrees: NDArray[np.float32], wavelength_nm: float, mds: MoonDatas, apply_correction: bool = True):
    a_l = np.exp(_ln_moon_disk_reflectance(absolute_mpa_degrees, wavelength_nm, mds))
    if apply_correction:
        mr_correction_factor = _get_correction_factor(
            wavelength_nm,
            np.radians(mds.ampa)
        )
        a_l = a_l * mr_correction_factor
    return a_l


def get_interpolated_reflectance(absolute_mpa_degrees: NDArray[np.float32], wavelength_nm: float, mds: MoonDatas, apply_correction: bool = True, adjust_apollo: bool = True):
    a_l = _interpolated_moon_disk_reflectance(absolute_mpa_degrees, wavelength_nm, mds, adjust_apollo)
    if apply_correction:
        mr_correction_factor = _get_correction_factor(
            wavelength_nm,
            np.radians(mds.ampa)
        )
        a_l = a_l * mr_correction_factor
    return a_l
