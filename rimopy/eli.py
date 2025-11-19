"""
Extraterrestrial Lunar Irradiance (ELI)

This module provides the main high-level interface for calculating the modeled
**Extraterrestrial Lunar Irradiance (ELI)**, following Román et al. (2020).

It allows users to compute the top-of-atmosphere lunar irradiance for any wavelength
and observation geometry, either from precomputed Moon data, from geographic coordinates
on Earth, or from custom observer kernels.

Exports
-------
Classes
    ELISettings
        Configuration controlling correction factors, Apollo adjustment, and output units.

Functions
    get_irradiance
        Main entry point to compute the modeled extraterrestrial lunar irradiance.
"""

from dataclasses import dataclass
from typing import List, Union, Iterable, overload

import numpy as np
from numpy.typing import NDArray

from . import esi, elref
from .types import MoonDatas, MissingRCFBehavior, EarthPoint
from .geometry import resolve_mds


@dataclass(frozen=True)
class ELISettings:
    """
    Settings of the methodology for calculating the Extraterrestrial Lunar Irradiance (ELI).

    Attributes
    ----------
    apply_correction : bool, default False
        If True, multiply the modeled irradiance by the RIMO Correction Factor (RCF) for
        each CIMEL wavelength (empirical correction from Román et al. 2020). If False,
        return the raw RIMO irradiance (no correction).
    adjust_apollo : bool, default True
        If True, adjust the ROLO model reflectance using Apollo spectra (derived from
        Apollo 16 lunar samples).
    per_nm : bool, default True
        If True, output ELI in W·m⁻²·nm⁻¹. Otherwise, output in W·m⁻².
    missing_rcf : MissingRCFBehavior, default MissingRCFBehavior.ERROR
        Behavior when at least one requested wavelength has no RCF available and
        `apply_correction` is True.
    """

    apply_correction: bool = False
    adjust_apollo: bool = True
    per_nm: bool = True
    missing_rcf: MissingRCFBehavior = MissingRCFBehavior.ERROR


def _get_esi(
    esi_calc: esi.ESICalculator,
    wavelengths_nm: Iterable[float],
    eli_settings: ELISettings,
) -> NDArray[np.float64]:
    """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
    Returns the data in Wm⁻²

    Parameters
    ----------
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    wavelengths_nm : iterable of float
        Wavelengths (in nanometers) of which the extraterrestrial solar irradiance will be
        obtained
    eli_settings : ELISettings
        Configuration of the ELI calculation method.

    Returns
    -------
    array of float
        The expected extraterrestrial solar irradiance in Wm⁻² or Wm⁻²/nm
    """
    return esi_calc.get_esi(wavelengths_nm, eli_settings.per_nm)


def _calculate_eli(
    wavelengths_nm: Iterable[float],
    mds: MoonDatas,
    esi_calc: esi.ESICalculator,
    eli_settings: ELISettings,
) -> NDArray[np.float64]:
    """Calculation of Extraterrestrial Lunar Irradiance following Eq 3 in Roman et al., 2020

    Simulates a lunar observation for a wavelength for any observer/solar selenographic
    latitude and longitude.

    Parameters
    ----------
    wavelength_nm : iterable of float
        Wavelengths (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated.
    mds : MoonDatas
        Moon data needed to calculate Moon's irradiance
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    eli_settings : ELISettings
        Configuration of the ELI calculation method.

    Returns
    -------
    array of float
        The extraterrestrial lunar irradiance calculated.
        One array per amount of moon geometry. Then, each inner array has the
        amount of values as the amount of wavelengths.
    """
    a_l = elref.get_reflectance(
        wavelengths_nm,
        mds=mds,
        apply_correction=eli_settings.apply_correction,
        missing_rcf=eli_settings.missing_rcf,
        adjust_apollo=eli_settings.adjust_apollo,
    )
    solid_angle_moon: float = 6.4177e-05
    omega = solid_angle_moon
    esk = _get_esi(esi_calc, wavelengths_nm, eli_settings)
    dsm = mds.dsm
    dom = mds.dom
    distance_earth_moon_km: int = 384400
    lunar_irr = (
        ((a_l * omega * esk) / np.pi).T
        * ((1 / dsm) ** 2)
        * (distance_earth_moon_km / dom) ** 2
    )
    return lunar_irr.T


@overload
def get_irradiance(
    wavelengths_nm: Iterable[float],
    *,
    mds: MoonDatas,
    esi_calc: esi.ESICalculator = None,
    eli_settings: ELISettings = None,
) -> NDArray[np.float64]:
    ...


@overload
def get_irradiance(
    wavelengths_nm: Iterable[float],
    *,
    earth_data: EarthPoint,
    kernels_path: str,
    esi_calc: esi.ESICalculator = None,
    eli_settings: ELISettings = None,
) -> NDArray[np.float64]:
    ...


@overload
def get_irradiance(
    wavelengths_nm: Iterable[float],
    *,
    utc_times: Union[str, List[str]],
    kernels_path: str,
    extra_kernels: List[str],
    extra_kernels_path: str,
    observer_name: str,
    esi_calc: esi.ESICalculator = None,
    eli_settings: ELISettings = None,
) -> NDArray[np.float64]:
    ...


def get_irradiance(
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
    esi_calc: esi.ESICalculator = None,
    eli_settings: ELISettings = None,
) -> NDArray[np.float64]:
    """
    Compute the Extraterrestrial Lunar Irradiance (ELI) following Eq. 3 in Román et al. (2020).

    This function calculates the modeled lunar irradiance at the top of the atmosphere for
    a given set of wavelengths and observation geometry.  The geometry can be provided in
    one of three equivalent ways:

      - **A)** `mds`: precomputed lunar geometry (`MoonDatas`).
      - **B)** `earth_data` + `kernels_path`: geographic coordinates and time(s) on Earth.
      - **C)** `utc_times`, `kernels_path`, `extra_kernels`, `extra_kernels_path`, `observer_name`:
        observer geometry defined by extra SPICE kernels.

    Parameters
    ----------
    wavelengths_nm : iterable of float
        Wavelengths in nanometers at which to compute the extraterrestrial lunar irradiance.
    mds : MoonDatas, optional
        Precomputed Moon geometry and distances (source A).
    earth_data : EarthPoint, optional
        Geographic location and times of the observation (source B).
    kernels_path : str, optional
        Directory containing the required SPICE kernels.
    utc_times : str or list of str, optional
        UTC datetimes of the observation(s) (source C).
    extra_kernels : list of str, optional
        Filenames of additional SPICE kernels describing the observer body (source C).
    extra_kernels_path : str, optional
        Directory containing the extra kernels (source C).
    observer_name : str, optional
        Name of the observer body as defined in the extra kernels (source C).
    esi_calc : esi.ESICalculator, optional
        Calculator for the extraterrestrial solar irradiance.
        Defaults to a linearly interpolated Wehrli-based implementation.
    eli_settings : ELISettings, optional
        Configuration controlling correction factors, Apollo adjustment, and output units.

    Returns
    -------
    ndarray of float
        Modeled extraterrestrial lunar irradiance in W·m⁻² (or W·m⁻²·nm⁻¹ if `per_nm=True` in
        `ELISettings`).  The output has shape ``(N_geometries, N_wavelengths)``; if there is only
        one geometry, the first dimension is squeezed and a one-dimensional array is returned.
    """
    if eli_settings is None:
        eli_settings = ELISettings()
    if esi_calc is None:
        esi_calc = esi.ESICalculatorWehrli()
    mds = resolve_mds(
        mds,
        earth_data,
        kernels_path,
        utc_times,
        extra_kernels,
        extra_kernels_path,
        observer_name,
    )
    irradiances = _calculate_eli(wavelengths_nm, mds, esi_calc, eli_settings)
    if len(irradiances) == 1:
        return irradiances[0]
    return irradiances
