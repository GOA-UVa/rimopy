"""
Extraterrestrial Lunar Reflectance (ELRef)

This module provides the functions to compute the modeled **Extraterrestrial Lunar Reflectance (ELRef)**,
following Eq. 2 in Román et al. (2020).  It evaluates the disk-integrated lunar reflectance at
given wavelengths and geometries, optionally applying RIMO correction factors and Apollo-based
spectral adjustments.

Exports
-------
Functions
    get_reflectance
        Main entry point to compute the modeled lunar disk reflectance.
"""

from typing import Iterable, Union, List, overload

import numpy as np
from numpy.typing import NDArray

from . import coefficients as coeffs
from .types import MoonDatas, MissingRCFBehavior, EarthPoint
from .geometry import resolve_mds
from . import correction_factor as corr_f


def _summatory_a(
    wavelengths_nm: NDArray[np.float64], gr: NDArray[np.float64]
) -> NDArray[np.float64]:
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
    wavelengths_nm: NDArray[np.float64], phi: NDArray[np.float64]
) -> NDArray[np.float64]:
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


def _ln_moon_disk_reflectance(
    wavelengths_nm: NDArray[np.float64],
    mds: MoonDatas,
) -> NDArray[np.float64]:
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
    wavelengths_nm: NDArray[np.float64],
    mds: MoonDatas,
    adjust_apollo: bool,
) -> NDArray[np.float64]:
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
    return np.array([np.interp(wavelengths_nm, x_values, yval) for yval in y_values])


@overload
def get_reflectance(
    wavelengths_nm: Iterable[float],
    *,
    mds: MoonDatas,
    apply_correction: bool = False,
    missing_rcf: MissingRCFBehavior = MissingRCFBehavior.ERROR,
    adjust_apollo: bool = True,
) -> NDArray[np.float64]:
    ...


@overload
def get_reflectance(
    wavelengths_nm: Iterable[float],
    *,
    earth_data: EarthPoint,
    kernels_path: str,
    apply_correction: bool = False,
    missing_rcf: MissingRCFBehavior = MissingRCFBehavior.ERROR,
    adjust_apollo: bool = True,
) -> NDArray[np.float64]:
    ...


@overload
def get_reflectance(
    wavelengths_nm: Iterable[float],
    *,
    utc_times: Union[str, List[str]],
    kernels_path: str,
    extra_kernels: List[str],
    extra_kernels_path: str,
    observer_name: str,
    apply_correction: bool = False,
    missing_rcf: MissingRCFBehavior = MissingRCFBehavior.ERROR,
    adjust_apollo: bool = True,
) -> NDArray[np.float64]:
    ...


def get_reflectance(
    wavelengths_nm: Iterable[float],
    *,
    # source A: directly moon datas
    mds: MoonDatas = None,
    # source B: earth point + kernels path
    earth_data: EarthPoint = None,
    # shared B & C
    kernels_path: str = None,
    # source C: extra kernels
    utc_times: Union[str, List[str]] = None,
    extra_kernels: List[str] = None,
    extra_kernels_path: str = None,
    observer_name: str = None,
    # common
    apply_correction: bool = False,
    missing_rcf: MissingRCFBehavior = MissingRCFBehavior.ERROR,
    adjust_apollo: bool = True,
) -> NDArray[np.float64]:
    """
    Compute the reflectance of the Moon's disk following Eq. 2 in Román et al. (2020).

    This function calculates the modeled disk-integrated lunar reflectance for one or more
    wavelengths. If the wavelength is not present in the ROLO coefficients, a linear interpolation
    (or extrapolation at the ends) is performed between the two nearest valid bands.

    Provide exactly one geometry source:

      - **A)** `mds`: precomputed Moon geometry data (`MoonDatas`).
      - **B)** `earth` + `kernels_path`: compute geometry from an EarthPoint and base SPICE kernels.
      - **C)** `utc_times` + `kernels_path` + `extra_kernels` + `extra_kernels_path` + `observer_name`:
        compute geometry using extra SPICE kernels for a custom observer body.

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths in nanometers for which to obtain the disk-integrated reflectance.
    mds : MoonDatas, optional
        Precomputed Moon geometry and distances.
    earth : EarthPoint, optional
        Geographic location and times of the observation (used with `kernels_path`).
    kernels_path : str, optional
        Directory containing the necessary SPICE kernels.
    utc_times : str or list of str, optional
        UTC datetimes of the observations (used with extra kernels).
    extra_kernels : list of str, optional
        List of extra kernel filenames defining the observer body.
    extra_kernels_path : str, optional
        Directory containing the extra kernels.
    observer_name : str, optional
        Name of the observer body as defined in the extra kernels.
    apply_correction : bool, default True
        If True, apply the RIMO Correction Factor (RCF) to the computed reflectance.
    missing_rcf : MissingRCFBehavior, default MissingRCFBehavior.ERROR
        Behavior when at least one requested wavelength lacks an RCF and
        `apply_correction` is True.
    adjust_apollo : bool, default True
        If True, adjust the modeled reflectance using Apollo spectral measurements.

    Returns
    -------
    ndarray of float
        The modeled lunar disk reflectance for the given geometries and wavelengths.
        The output has shape ``(N_geometries, N_wavelengths)``.
    """
    mds = resolve_mds(
        mds,
        earth_data,
        kernels_path,
        utc_times,
        extra_kernels,
        extra_kernels_path,
        observer_name,
    )
    wavelengths_nm = np.array(wavelengths_nm)
    a_l = _interpolated_moon_disk_reflectance(wavelengths_nm, mds, adjust_apollo)
    if apply_correction:
        correction_factor = corr_f.get_correction_factor(
            wavelengths_nm,
            np.radians(mds.mpa),
            missing_rcf,
        )
        a_l = a_l * correction_factor
    return a_l
