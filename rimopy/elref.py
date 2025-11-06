"""ELRef Extraterrestrial Lunar Reflectance

This module allows for the calculation of the extraterrestrial lunar reflectance.
"""

from typing import Iterable

import numpy as np
from numpy.typing import NDArray

from . import coefficients as coeffs
from .types import MoonDatas
from . import correction_factor as corr_f


def _summatory_a(
    wavelengths_nm: NDArray[np.float32], gr: NDArray[np.float32]
) -> NDArray[np.float32]:
    """The first summatory of Eq. 2 in Roman et al., 2020

    Parameters
    ----------
    wavelength_nm : array of float
        Wavelengths in nanometers from which the moon's disk reflectance is being calculated
    gr : array of float
        Absolute value of MPA in radians

    Returns
    -------
    array of float
        Result of the computation of the first summatory. One array per amount of `gr`.
        Then, each inner array has the amount of values as the amount of wavelengths.
    """
    ac = coeffs.get_coefficients_a(wavelengths_nm)
    gr = np.array([gr]).T
    sa = ac[0] + ac[1] * gr + ac[2] * gr**2 + ac[3] * gr**3
    return sa


def _summatory_b(
    wavelengths_nm: NDArray[np.float32], phi: NDArray[np.float32]
) -> NDArray[np.float32]:
    """The second summatory of Eq. 2 in Roman et al., 2020, without the erratum

    Parameters
    ----------
    wavelengths_nm : array of float
        Wavelengths from which the moon's disk reflectance is being calculated
    phi : array of float
        Selenographic longitude of the Sun (in radians)

    Returns
    -------
    array of float
        Result of the computation of the second summatory. One array per amount of `phi`.
        Then, each inner array has the amount of values as the amount of wavelengths.
    """
    bc = coeffs.get_coefficients_b(wavelengths_nm)
    phi = np.array([phi]).T
    sb = bc[0] * phi + bc[1] * phi**3 + bc[2] * phi**5
    return sb


def _get_correction_factor(
    wavelengths_nm: Iterable[float], mpa: NDArray[np.float32]
) -> NDArray[np.float32]:
    """Calculation of RIMO correction factor (RCF) following Eq 9 in Roman et al., 2020

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated
    mpa : array of float
        Absolute Moon phase angle (in radians)

    Returns
    -------
    array of float
        The calculated RCF. One array per amount of `mpa`.
        Then, each inner array has the amount of values as the amount of wavelengths.
    """
    params = corr_f.get_correction_params(wavelengths_nm)
    mpa = np.array([mpa]).T
    rcf = params.a_coeff + params.b_coeff * mpa + params.c_coeff * mpa**2
    return rcf


def _ln_moon_disk_reflectance(
    wavelengths_nm: NDArray[np.float32],
    mds: MoonDatas,
) -> NDArray[np.float32]:
    """The calculation of the ln of the reflectance of the Moon's disk, following Eq.2 in
    Roman et al., 2020

    If the wavelength has no associated ROLO coefficients, it uses some linearly interpolated
    ones.

    Parameters
    ----------
    wavelength_nm : array of float
        Wavelengths in nanometers from which one wants to obtain the MDRs.
    mds : MoonDatas
        Moon data needed to calculate Moon's irradiance

    Returns
    -------
    array of float
        The ln of the reflectance of the Moon's disk for the inputed data.
        One array per amount of moon geometry. Then, each inner array has the
        amount of values as the amount of wavelengths.
    """
    ampa = mds.ampa
    gr_value = np.radians(ampa)
    phi = mds.lonsun
    c_coeffs = coeffs.get_coefficients_c()
    d_coeffs = coeffs.get_coefficients_d(wavelengths_nm)
    p_coeffs = coeffs.get_coefficients_p()
    sum_a = _summatory_a(wavelengths_nm, gr_value)
    sum_b = _summatory_b(wavelengths_nm, phi)
    gd_value = np.array([ampa]).T
    d1_value = d_coeffs[0] * np.exp(-gd_value / p_coeffs[0])
    d2_value = d_coeffs[1] * np.exp(-gd_value / p_coeffs[1])
    d3_value = d_coeffs[2] * np.cos((gd_value - p_coeffs[2]) / p_coeffs[3])
    phi = np.array([phi]).T
    l_theta = np.array([mds.latobs]).T
    l_phi = np.array([mds.lonobs]).T
    result = (
        sum_a
        + sum_b
        + c_coeffs[0] * l_phi
        + c_coeffs[1] * l_theta
        + c_coeffs[2] * phi * l_phi
        + c_coeffs[3] * phi * l_theta
        + d1_value
        + d2_value
        + d3_value
    )
    return result


def _neighbors_set_linear_exact(query: Iterable[float], reference: Iterable[float]):
    """
    Find nearest neighbors of query values within a sorted reference sequence.

    For each element in `query`, the function selects:
      - The value itself, if it exists in `reference`.
      - Otherwise, the nearest left and right neighbors in `reference`.
      - If the query value is out of the reference range, the two closest
        boundary values are included.

    Parameters
    ----------
    query : Iterable of float
        Values for which to find nearest neighbors.
    reference : Iterable of float
        Sorted sequence of reference values (must be in ascending order).

    Returns
    -------
    neighbors : set of float
        Unique set of reference values that are either exact matches or
        nearest neighbors of the query values.
    """
    n_ref = len(reference)
    if n_ref == 0:
        return set()
    if n_ref == 1:
        return set(reference)
    query = sorted(query)
    j = 0
    out = set()
    first_two = (reference[0], reference[1])
    last_two = (reference[-2], reference[-1])
    for qv in query:
        # advance j until reference[j] >= qv (or j == n_ref)
        while j < n_ref and reference[j] < qv:
            j += 1
        if j < n_ref and reference[j] == qv:
            # exact hit: include only qv
            out.add(reference[j])
        elif j == 0:
            # qv < reference[0]
            out.add(first_two[0])
            out.add(first_two[1])
        elif j == n_ref:
            # qv > reference[-1]
            out.add(last_two[0])
            out.add(last_two[1])
        else:
            # interior miss: include neighbors
            out.add(reference[j - 1])
            out.add(reference[j])
    return out


def _interpolated_moon_disk_reflectance(
    wavelengths_nm: NDArray[np.float32],
    mds: MoonDatas,
    adjust_apollo: bool,
) -> NDArray[np.float32]:
    """The calculation of the reflectance of the Moon's disk, following Eq.2 in Roman et al., 2020

    If the wavelength is not present in the ROLO coefficients, it calculates the linear
    interpolation between the previous and the next one, or the extrapolation with the two
    nearest ones in case that it's on an extreme.

    Parameters
    ----------
    wavelengths_nm : array of float
        Wavelengths in nanometers from which one wants to obtain the MDR.
    mds : MoonDatas
        Moon data needed to calculate Moon's irradiance
    adjust_apollo : bool
        If True, the calculated reflectance will be adjusted to the Apollo spectra.

    Returns
    -------
    array of float
        The ln of the reflectance of the Moon's disk for the inputed data.
        One array per amount of moon geometry. Then, each inner array has the
        amount of values as the amount of wavelengths.
    """
    wvlens = coeffs.get_wavelengths()
    if adjust_apollo:
        apollo_coeffs = coeffs.get_apollo_coefficients()
    else:
        apollo_coeffs = np.ones(shape=len(wvlens))
    wavelengths_nm = np.where(wavelengths_nm < wvlens[0], wvlens[0], wavelengths_nm)
    wavelengths_nm = np.where(wavelengths_nm > wvlens[-1], wvlens[-1], wavelengths_nm)
    x_values = _neighbors_set_linear_exact(wavelengths_nm, wvlens)
    x_values = np.array(sorted(x_values))
    ap_indices = np.where(np.isin(wvlens, x_values))[0]
    y_values = (
        np.exp(_ln_moon_disk_reflectance(x_values, mds)) * apollo_coeffs[ap_indices]
    )
    return np.array([np.interp(wavelengths_nm, x_values, yval) for yval in y_values]).T


def get_reflectance_interpolating_coefficients(
    wavelengths_nm: Iterable[float],
    mds: MoonDatas,
    apply_correction: bool = True,
):
    """The calculation of the reflectance of the Moon's disk, following Eq.2 in Roman et al., 2020

    If the wavelength has no associated ROLO coefficients, it uses some linearly interpolated
    ones.

    Parameters
    ----------
    wavelength_nm : iterable of float
        Wavelengths in nanometers from which one wants to obtain the MDRs.
    mds : MoonDatas
        Moon data needed to calculate Moon's irradiance
    apply_correction: bool
        If True the RIMO Correction Factor will be calculated and applied to the obtained
        reflectance.

    Returns
    -------
    array of float
        The ln of the reflectance of the Moon's disk for the inputed data.
        One array per amount of moon geometry. Then, each inner array has the
        amount of values as the amount of wavelengths.
    """
    wavelengths_nm = np.array(wavelengths_nm)
    a_l = np.exp(_ln_moon_disk_reflectance(wavelengths_nm, mds))
    if apply_correction:
        mr_correction_factor = _get_correction_factor(
            wavelengths_nm, np.radians(mds.ampa)
        ).T
        a_l = a_l * mr_correction_factor
    return a_l


def get_interpolated_reflectance(
    wavelengths_nm: Iterable[float],
    mds: MoonDatas,
    apply_correction: bool = True,
    adjust_apollo: bool = True,
):
    """The calculation of the reflectance of the Moon's disk, following Eq.2 in Roman et al., 2020

    If the wavelength is not present in the ROLO coefficients, it calculates the linear
    interpolation between the previous and the next one, or the extrapolation with the two
    nearest ones in case that it's on an extreme.

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths in nanometers from which one wants to obtain the MDR.
    mds : MoonDatas
        Moon data needed to calculate Moon's irradiance
    apply_correction: bool
        If True the RIMO Correction Factor will be calculated and applied to the obtained
        reflectance.
    adjust_apollo : bool
        If True, the calculated reflectance will be adjusted to the Apollo spectra.

    Returns
    -------
    array of float
        The ln of the reflectance of the Moon's disk for the inputed data.
        One array per amount of moon geometry. Then, each inner array has the
        amount of values as the amount of wavelengths.
    """
    wavelengths_nm = np.array(wavelengths_nm)
    a_l = _interpolated_moon_disk_reflectance(wavelengths_nm, mds, adjust_apollo)
    if apply_correction:
        mr_correction_factor = _get_correction_factor(
            wavelengths_nm, np.radians(mds.ampa)
        ).T
        a_l = a_l * mr_correction_factor
    return a_l
