"""Correction Factor

This module contains the coefficients of the RIMO correction factor (RCF).

This value corrects the simulated extra terrestrial lunar irradiance in order
to use it with photometers.

See "Roman et al., 2020: Correction of a lunar-irradiance model for aerosol optical depth
retrieval and comparison with a star photometer" for mor information.

It exports the following classes:
    * CorrectionParams - DataClass that contains the estimated coefficients of the
        RCF for a wavelength.

It exports the following functions:

    * get_correction_params - returns the RCF coefficients estimated for a wavelength
"""

from dataclasses import dataclass
from typing import List, Iterable
import logging

from numpy.typing import NDArray
import numpy as np

from .types import MissingRCFBehavior


@dataclass
class CorrectionParams:
    """
    DataClass that contains the estimated coefficients of the RCF for a wavelength.

    Attributes
    ----------
    a_coeff : np.array of float
        RCF coefficient 'a'
    b_coeff : np.array of float
        RCF coefficient 'b'
    c_coeff : np.array of float
        RCF coefficient 'c'
    """

    a_coeff: NDArray[np.float32]
    b_coeff: NDArray[np.float32]
    c_coeff: NDArray[np.float32]


def _get_corrected_wavelengths() -> List[float]:
    """Gets all wavelengths (in nanometers) presented in the RCF model

    Returns
    -------
    list of float
        A list of floats that are the wavelengths in nanometers, in order
    """
    return [340, 380, 440, 500, 675, 870, 935, 1020, 1640]


def _get_all_correction_params() -> List[List[float]]:
    """Gets all RCF coefficients

    Returns
    -------
    list of list of float
        A list containing multiple list of floats. Each sublist is the list of coefficients
        for a wavelength
    """
    return [
        [1.186, -2.35 * 10**-2, 1.92 * 10**-1],
        [1.082, -4.17 * 10**-3, 7.10 * 10**-2],
        [1.062, -5.35 * 10**-4, 1.14 * 10**-2],
        [1.078, -8.93 * 10**-4, 1.11 * 10**-2],
        [1.092, -4.50 * 10**-4, 1.38 * 10**-2],
        [1.075, -2.05 * 10**-3, 1.37 * 10**-2],
        [1.071, -2.41 * 10**-3, 1.36 * 10**-2],
        [1.035, 5.55 * 10**-3, 2.79 * 10**-2],
        [1.047, -1.25 * 10**-3, 2.26 * 10**-2],
    ]


def _get_all_as() -> List[float]:
    """Gets all 'a' RCF coefficients

    Returns
    -------
    list of float
        A list containing all 'a' coefficients in wavelength order
    """
    return list(map(lambda x: x[0], _get_all_correction_params()))


def _get_all_bs() -> List[float]:
    """Gets all 'b' RCF coefficients

    Returns
    -------
    list of float
        A list containing all 'b' coefficients in wavelength order
    """
    return list(map(lambda x: x[1], _get_all_correction_params()))


def _get_all_cs() -> List[float]:
    """Gets all 'c' RCF coefficients

    Returns
    -------
    list of float
        A list containing all 'c' coefficients in wavelength order
    """
    return list(map(lambda x: x[2], _get_all_correction_params()))


def _get_correction_params_fill_ones(
    wavelengths_nm: Iterable[float], atol_nm: float
) -> "CorrectionParams":
    """Obtain the RCF params, and mock values for invalid wavelengths

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths (in nanometers) of which one wants to obtain the RCF params.
        If there's no RCF for the wavelength, it will return (1, 0, 0), which will
        render a RCF = 1.
    atol_nm: float
        Absolute tolerance for matching wavelengths to the supported set.

    Returns
    -------
    'CorrectionParams'
        RIMO RCF correction params
    """
    x_values = _get_corrected_wavelengths()
    wl = np.array(wavelengths_nm)
    all_as = _get_all_as()
    all_bs = _get_all_bs()
    all_cs = _get_all_cs()
    a_coeff = np.ones_like(wl, dtype=float)
    b_coeff = np.zeros_like(wl, dtype=float)
    c_coeff = np.zeros_like(wl, dtype=float)
    # Fill real coefficients ONLY for wavelengths that have an RCF defined.
    # Treat "has RCF" as "equals one of x_values within a tiny tolerance".
    for i, w in enumerate(wl):
        idx = np.where(np.isclose(x_values, w, rtol=0.0, atol=atol_nm))[0]
        if idx.size:  # wavelength is supported → copy its coefficients
            j = int(idx[0])
            a_coeff[i] = all_as[j]
            b_coeff[i] = all_bs[j]
            c_coeff[i] = all_cs[j]
    return CorrectionParams(a_coeff, b_coeff, c_coeff)


def _get_correction_params(
    wavelengths_nm: Iterable[float],
    missing_rcf: MissingRCFBehavior,
    atol_nm: float = 0.1,
) -> "CorrectionParams":
    """Gets the RCF correction parameters for a specific wavelength in nanometers

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths (in nanometers) of which one wants to estimate the RCF params
    missing_rcf: MissingRCFBehavior
        Behavior when at least one requested wavelength has no RCF available.
    atol_nm : float, default 0.1
        Absolute tolerance for matching wavelengths to the supported set.

    Returns
    -------
    'CorrectionParams'
        Arrays of (a, b, c) aligned with the input order.
    """
    wavelengths_nm = np.array(list(wavelengths_nm))
    supported = np.array(_get_corrected_wavelengths())
    mask = np.any(
        np.isclose(wavelengths_nm[:, None], supported[None, :], rtol=0.0, atol=atol_nm),
        axis=1,
    )
    if not np.all(mask):
        missing_list: List[float] = [float(w) for w in wavelengths_nm[~mask]]
        msg = f"RCF not available for the wavelengths: {missing_list}"
        if missing_rcf is MissingRCFBehavior.ERROR:
            raise ValueError(msg)
        elif missing_rcf is MissingRCFBehavior.WARN:
            logging.warning(msg)
    pars = _get_correction_params_fill_ones(wavelengths_nm, atol_nm)
    return pars


def _calc_correction_factor(
    wavelengths_nm: Iterable[float],
    mpa: NDArray[np.float32],
    missing_rcf: MissingRCFBehavior,
    atol_nm: float = 0.1,
) -> NDArray[np.float32]:
    """Calculation of RIMO correction factor (RCF) following Eq 9 in Roman et al., 2020

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths (in nanometers) of which the RCFs will be calculated
    mpa : array of float
        Moon phase angle (in radians)
    missing_rcf: MissingRCFBehavior
        Behavior when at least one requested wavelength has no RCF available and
        `apply_correction` is True.
    Returns
    -------
    array of float
        The calculated RCF. One array per amount of `mpa`.
        Then, each inner array has the amount of values as the amount of wavelengths.

    """
    params = _get_correction_params(wavelengths_nm, missing_rcf, atol_nm)
    mpa = np.array([mpa]).T
    rcf = params.a_coeff + params.b_coeff * mpa + params.c_coeff * mpa**2
    return rcf


def _nearest_rcfs(wavelengths_nm: Iterable[float], mpa: NDArray[np.float32]):
    """Obtain the RCF of the closest wavelength if the wavelength has no RCF available

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths (in nanometers) of which one wants to obtain the RCF
    mpa : array of float
        Moon phase angle (in radians)

    Returns
    -------
    array of float
        The estimated RCF. One array per amount of `mpa`.
        Then, each inner array has the amount of values as the amount of wavelengths.
    """
    x_values = np.array(_get_corrected_wavelengths())
    y_values = _calc_correction_factor(x_values, mpa, MissingRCFBehavior.IGNORE)
    wavelengths_nm = np.asarray(wavelengths_nm, dtype=float)
    wavelengths_nm = np.clip(wavelengths_nm, x_values[0], x_values[-1])
    wavelengths_nm = np.where(wavelengths_nm < x_values[0], x_values[0], wavelengths_nm)
    wavelengths_nm = np.where(
        wavelengths_nm > x_values[-1], x_values[-1], wavelengths_nm
    )
    diff = np.abs(x_values[:, None] - wavelengths_nm[None, :])
    nearest_idx = diff.argmin(axis=0)
    return y_values[:, nearest_idx]


def get_correction_factor(
    wavelengths_nm: Iterable[float],
    mpa: NDArray[np.float32],
    missing_rcf: MissingRCFBehavior,
) -> NDArray[np.float32]:
    """Calculation of RIMO correction factor (RCF) following Eq 9 in Roman et al., 2020

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated
    mpa : array of float
        Moon phase angle (in radians)
    missing_rcf: MissingRCFBehavior
        Behavior when at least one requested wavelength has no RCF available and
        `apply_correction` is True.
    Returns
    -------
    array of float
        The calculated RCF. One array per amount of `mpa`.
        Then, each inner array has the amount of values as the amount of wavelengths.
    """
    if missing_rcf == MissingRCFBehavior.NEAREST:
        rcf = _nearest_rcfs(wavelengths_nm, mpa)
    else:
        rcf = _calc_correction_factor(wavelengths_nm, mpa, missing_rcf)
    return rcf
