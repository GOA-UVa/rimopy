"""ELI Extraterrestrial Lunar Irradiance

This module is the main module, as it allowes the user to calculate the Extraterrestrial Lunar
Irradiance at a concrete wavelength, at an absolute Moon phase angle, and giving selenographic
parameters.

It exports the following classes:
    * ELISettings - Settings that will modify the methodology of calculating the ELI
    * EarthPoint - Data of the point on Earth surface of which the ELI will be calculated.

It exports the following functions:

    * get_eli_bypass - returns the expected extraterrestrial lunar irradiation of a wavelength for
        any observer/solar selenographic coordinates, in Wm⁻².
    * get_eli - returns the expected extraterrestrial lunar irradiation of a wavelength in any
        geographic coordinates, in Wm⁻².
    * get_eli_from_extra_kernels - returns the expected extraterrestrial lunar irradiation of
        a wavelength in any geographic coordinates, in Wm⁻², using data from extra kernels for the
        observer body.
    * get_eli_bypass_per_nm - returns the expected extraterrestrial lunar irradiation of a
        wavelength for any observer/solar selenographic coordinates, in Wm⁻²/nm.
    * get_eli_per_nm - returns the expected extraterrestrial lunar irradiation of a wavelength in
        any geographic coordinates, in Wm⁻²/nm.
    * get_eli_per_nm_from_extra_kernels - returns the expected extraterrestrial lunar irradiation of
        a wavelength in any geographic coordinates, in Wm⁻²/nm, using data from extra kernels for the
        observer body.
"""

from dataclasses import dataclass
import math
from typing import List, Union

from . import spice_iface, esi, elref
from .types import MoonData

@dataclass(frozen=True)
class ELISettings:
    """
    Settings that will modify the methodology of calculating the ELI

    Attributes
    ----------
    apply_correction : bool
        If True the result will have been multiplied by the RCF (Rimo Correction Factor),
        which corrects the data for the photometers' calibration.
        Otherwise it won't.
    interpolate_rolo_coefficients : bool
        If True the reflectance will be calculated linearly interpolating the ROLO coefficients.
        Otherwise it will be calculated interpolating the surrounding reflectances, calculated
        with empirical coefficients.
    adjust_apollo : bool
        If True the ROLO model reflectance will be adjusted using Apollo spectra, in case it's
        calculated interpolating surrounding reflectances. The Apollo spectra is the spectra
        generated with the Moon samples from Apollo 16th mission.
    """
    apply_correction: bool = False
    interpolate_rolo_coefficients: bool = False
    adjust_apollo: bool = True

@dataclass
class EarthPoint:
    """
    Data of the point on Earth surface of which the ELI will be calculated.

    Attributes
    ----------
    lat : float
        Geographic latitude (in degrees) of the location.
    lon : float
        Geographic longitude (in degrees) of the location.
    utc_times : list of str | str
        Time/s at which the ELI will be calculated, in a valid UTC DateTime format.
    altitude : float
        Altitude over the sea level in meters. Default = 0.
    """
    __slots__ = ['lat', 'lon', 'utc_times', 'altitude']
    def __init__(self, lat: float, lon: float, utc_times: Union[List[str], str],
                 altitude: float = 0):
        """
        Parameters
        ----------
        lat : float
            Geographic latitude (in degrees) of the location.
        lon : float
            Geographic longitude (in degrees) of the location.
        utc_times : list of str | str
            Time/s at which the ELI will be calculated, in a valid UTC DateTime format.
        altitude : float
            Altitude over the sea level in meters. Default = 0.
        """
        self.lat = lat
        self.lon = lon
        self.altitude = altitude
        if isinstance(utc_times, list):
            self.utc_times = utc_times
        else:
            self.utc_times = [utc_times]

    def set_utc_times(self, utc_times: Union[List[str], str]):
        """
        Modifies the utc_times attribute

        Parameters
        ----------
        utc_times : list of str | str
            Time/s at which the ELI will be calculated, in a valid UTC DateTime format.
        """
        if isinstance(utc_times, list):
            self.utc_times = utc_times
        else:
            self.utc_times = [utc_times]


def _get_esi(esi_calc: esi.ESICalculator, wavelength_nm: float) -> float:
    """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
    Returns the data in Wm⁻²

    Parameters
    ----------
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial solar irradiance will be
        obtained

    Returns
    -------
    float
        The expected extraterrestrial solar irradiance in W/sm
    """
    return esi_calc.get_esi(wavelength_nm)

def _get_esi_per_nm(esi_calc: esi.ESICalculator, wavelength_nm: float) -> float:
    """Gets the expected extraterrestrial solar irradiance at a concrete wavelength
    Returns the data in Wm⁻²/nm

    Parameters
    ----------
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial solar irradiance will be
        obtained

    Returns
    -------
    float
        The expected extraterrestrial solar irradiance in W/sm
        """
    return esi_calc.get_esi_per_nm(wavelength_nm)

def _calculate_eli(wavelength_nm: float, moon_data: MoonData, esi_calc: esi.ESICalculator,
                   eli_settings: ELISettings, per_nm: bool = False) -> float:
    """Calculation of Extraterrestrial Lunar Irradiance following Eq 3 in Roman et al., 2020

    Simulates a lunar observation for a wavelength for any observer/solar selenographic
    latitude and longitude.

    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated.
    moon_data : MoonData
        Moon data needed to calculate Moon's irradiance
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    eli_settings : ELISettings
        Configuration of the ELI calculation method.
    per_nm : bool
        True if the user wants the ELI in Wm⁻²/nm, otherwise it will be in Wm⁻². Default is False.

    Returns
    -------
    float
        The extraterrestrial lunar irradiance calculated
    """
    if not eli_settings.interpolate_rolo_coefficients:
        a_l = elref.get_interpolated_reflectance(
            moon_data.absolute_mpa_degrees,
            wavelength_nm,
            moon_data,
            eli_settings.apply_correction,
            eli_settings.adjust_apollo
        )
    else:
        a_l = elref.get_reflectance_interpolating_coefficients(
            moon_data.absolute_mpa_degrees,
            wavelength_nm,
            moon_data,
            eli_settings.apply_correction,
        )

    solid_angle_moon: float = 6.4177e-05
    omega = solid_angle_moon
    if per_nm:
        esk = _get_esi_per_nm(esi_calc, wavelength_nm)
    else:
        esk = _get_esi(esi_calc, wavelength_nm)
    dsm = moon_data.distance_sun_moon
    dom = moon_data.distance_observer_moon
    distance_earth_moon_km: int = 384400

    lunar_irr = ((a_l * omega * esk) / math.pi) * ((1 / dsm) ** 2) * \
        (distance_earth_moon_km / dom) ** 2
    return lunar_irr

def get_eli_bypass(wavelength_nm: Union[float, List[float]], moon_data: MoonData,
                   esi_calc: esi.ESICalculator = None,
                   eli_settings: ELISettings = None) -> Union[float, List[float]]:
    """Calculation of Extraterrestrial Lunar Irradiance following Eq 3 in Roman et al., 2020

    Allow users to simulate lunar observation for any observer/solar selenographic
    latitude and longitude (thus bypassing the need for their computation from the
    position/time of the observer).

    Returns the data in Wm⁻²

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated.
    moon_data : MoonData
        Moon data needed to calculate Moon's irradiance
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    eli_settings : ELISettings
        Configuration of the ELI calculation method.

    Returns
    -------
    float | list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter
        "wavelength_nm" was a list.
    """
    if eli_settings is None:
        eli_settings = ELISettings()
    if esi_calc is None:
        esi_calc = esi.ESICalculator()
    if isinstance(wavelength_nm, list):
        elis = []
        for wlen in wavelength_nm:
            elis.append(_calculate_eli(wlen, moon_data, esi_calc, eli_settings))
        return elis
    return _calculate_eli(wavelength_nm, moon_data, esi_calc, eli_settings)

def get_eli_from_extra_kernels(wavelength_nm: Union[float, List[float]],
        utc_times: Union[str, List[str]], kernels_path: str, extra_kernels: List[str],
        extra_kernels_path: str, observer_name: str,
        esi_calc: esi.ESICalculator = esi.ESICalculator(),
        eli_settings: ELISettings = ELISettings()
        ) -> Union[float, List[float], List[List[float]]]:
    """Calculation of Extraterrestrial Lunar Irradiance from geographic coordinates

    Allow users to simulate lunar observations for any observer position around the Earth
    and at any time.

    It loads the observer body data from custom extra kernels instead of generating it from
    basic kernels.

    Returns the data in Wm⁻²

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated.
    utc_times: str | list of str
        Time/s at which the ELI will be calculated, in a valid UTC DateTime format.
    kernels_path : str
        Folder where the needed SPICE kernels are stored.
    extra_kernels: list of str
        Custom kernels from which the observer body will be loaded, instead of calculating it.
    extra_kernels_path: str
        Folder where the extra kernels are located.
    observer_name: str
        Name of the body of the observer that will be loaded from the extra kernels.
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    eli_settings : ELISettings
        Configuration of the ELI calculation method.

    Returns
    -------
    float | list of float | list of list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter
        "wavelength_nm" was a list, or if EarthPoint utc_times is a list. In case both are lists,
        it will be a list of lists of floats. The list will contain a list for every utc_time
        present, and each of those will contain all irradiances associated for each wavelength.
    """
    moon_data = spice_iface.get_moon_datas_from_extra_kernels(utc_times, kernels_path,
                                                              extra_kernels, extra_kernels_path,
                                                              observer_name)
    irradiances = []
    for moon_d in moon_data:
        irradiances.append(get_eli_bypass(wavelength_nm, moon_d, esi_calc, eli_settings))
    if len(irradiances) == 1:
        return irradiances[0]
    return irradiances

def get_eli(wavelength_nm: Union[float, List[float]], earth_data: EarthPoint, kernels_path: str,
            esi_calc: esi.ESICalculator = None,
            eli_settings: ELISettings = None,
            ) -> Union[float, List[float], List[List[float]]]:
    """Calculation of Extraterrestrial Lunar Irradiance from geographic coordinates

    Allow users to simulate lunar observations for any observer position around the Earth
    and at any time.

    Returns the data in Wm⁻²

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated.
    earth_data : EarthPoint
        Data of the point on Earth surface of which the ELI will be calculated.
    kernels_path : str
        Folder where the needed SPICE kernels are stored.
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    eli_settings : ELISettings
        Configuration of the ELI calculation method.

    Returns
    -------
    float | list of float | list of list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter
        "wavelength_nm" was a list, or if EarthPoint utc_times is a list. In case both are lists,
        it will be a list of lists of floats. The list will contain a list for every utc_time
        present, and each of those will contain all irradiances associated for each wavelength.
    """
    moon_data = spice_iface.get_moon_datas(earth_data.lat, earth_data.lon, earth_data.altitude,
                                           earth_data.utc_times, kernels_path)
    irradiances = []
    for moon_d in moon_data:
        irradiances.append(get_eli_bypass(wavelength_nm, moon_d, esi_calc, eli_settings))
    if len(irradiances) == 1:
        return irradiances[0]
    return irradiances

def get_eli_bypass_per_nm(wavelength_nm: Union[float, List[float]], moon_data: MoonData,
                          esi_calc: esi.ESICalculator = None,
                          eli_settings: ELISettings = None,
                          ) -> Union[float, List[float]]:
    """Calculation of Extraterrestrial Lunar Irradiance following Eq 3 in Roman et al., 2020

    Allow users to simulate lunar observation for any observer/solar selenographic
    latitude and longitude (thus bypassing the need for their computation from the
    position/time of the observer).

    Returns the data in Wm⁻²/nm

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated.
    moon_data : MoonData
        Moon data needed to calculate Moon's irradiance.
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    eli_settings : ELISettings
        Configuration of the ELI calculation method.

    Returns
    -------
    float | list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter
        "wavelength_nm" was a list.
    """
    if eli_settings is None:
        eli_settings = ELISettings()
    if esi_calc is None:
        esi_calc = esi.ESICalculator()
    if isinstance(wavelength_nm, list):
        elis = []
        for wlen in wavelength_nm:
            elis.append(_calculate_eli(wlen, moon_data, esi_calc, eli_settings, True))
        return elis
    return _calculate_eli(wavelength_nm, moon_data, esi_calc, eli_settings, True)

def get_eli_per_nm_from_extra_kernels(wavelength_nm: Union[float, List[float]],
        utc_times: Union[str, List[str]], kernels_path: str, extra_kernels: List[str],
        extra_kernels_path: str, observer_name: str,
        esi_calc: esi.ESICalculator = None,
        eli_settings: ELISettings = None,
        ) -> Union[float, List[float], List[List[float]]]:
    """Calculation of Extraterrestrial Lunar Irradiance from geographic coordinates

    Allow users to simulate lunar observations for any observer position around the Earth
    and at any time.

    It loads the observer body data from custom extra kernels instead of generating it from
    basic kernels.

    Returns the data in Wm⁻²/nm

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated.
    utc_times: str | list of str
        Time/s at which the ELI will be calculated, in a valid UTC DateTime format.
    kernels_path : str
        Folder where the needed SPICE kernels are stored.
    extra_kernels: list of str
        Custom kernels from which the observer body will be loaded, instead of calculating it.
    extra_kernels_path: str
        Folder where the extra kernels are located.
    observer_name: str
        Name of the body of the observer that will be loaded from the extra kernels.
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    eli_settings : ELISettings
        Configuration of the ELI calculation method.

    Returns
    -------
    float | list of float | list of list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter
        "wavelength_nm" was a list, or if EarthPoint utc_times is a list. In case both are lists,
        it will be a list of lists of floats. The list will contain a list for every utc_time
        present, and each of those will contain all irradiances associated for each wavelength.
    """
    moon_data = spice_iface.get_moon_datas_from_extra_kernels(utc_times, kernels_path,
                                                              extra_kernels, extra_kernels_path,
                                                              observer_name)
    irradiances = []
    for moon_d in moon_data:
        irradiances.append(get_eli_bypass_per_nm(wavelength_nm, moon_d, esi_calc, eli_settings))
    if len(irradiances) == 1:
        return irradiances[0]
    return irradiances

def get_eli_per_nm(wavelength_nm: Union[float, List[float]], earth_data: EarthPoint,
                   kernels_path: str, esi_calc: esi.ESICalculator = None,
                   eli_settings: ELISettings = None,
                   ) -> Union[float, List[float], List[List[float]]]:
    """Calculation of Extraterrestrial Lunar Irradiance from geographic coordinates

    Allow users to simulate lunar observations for any observer position around the Earth
    and at any time.

    Returns the data in Wm⁻²/nm

    Parameters
    ----------
    wavelength_nm : float | list of float
        Wavelength/s (in nanometers) of which the extraterrestrial lunar irradiance will be
        calculated.
    earth_data : EarthPoint
        Data of the point on Earth surface of which the ELI will be calculated.
    kernels_path : str
        Folder where the needed SPICE kernels are stored.
    esi_calc : esi.ESICalculator
        ESI Calculator that will be used in the calculation of the Extraterrestrial Solar
        Irradiance.
    eli_settings : ELISettings
        Configuration of the ELI calculation method.

    Returns
    -------
    float | list of float | list of list of float
        The extraterrestrial lunar irradiance/s calculated. It will be a list if parameter
        "wavelength_nm" was a list, or if EarthPoint utc_times is a list. In case both are lists,
        it will be a list of lists of floats. The list will contain a list for every utc_time
        present, and each of those will contain all irradiances associated for each wavelength.
    """
    moon_data = spice_iface.get_moon_datas(earth_data.lat, earth_data.lon,
                                               earth_data.altitude, earth_data.utc_times,
                                               kernels_path)
    irradiances = []
    for moon_d in moon_data:
        irradiances.append(get_eli_bypass_per_nm(wavelength_nm, moon_d, esi_calc, eli_settings))
    if len(irradiances) == 1:
        return irradiances[0]
    return irradiances
