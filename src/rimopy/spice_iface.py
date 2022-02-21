"""SPICE iface

Interface with NASA's SPICE toolbox

It exports the following functions:

    * get_moon_datas - Calculates needed MoonData from SPICE toolbox
"""

from dataclasses import dataclass
import os
import math
from typing import List
import numpy as np
import spiceypy as spice
from .types import MoonData

CUSTOM_KERNEL_NAME = "custom.bsp"
EARTH_ID_CODE = 399

def _calculate_states(ets: np.ndarray, pos_iau_earth: np.ndarray, delta_t: float,
                      source_frame: str, target_frame: str) -> np.ndarray:
    """
    Returns a ndarray containing the states of a point referencing the target frame.

    The states array is a time-ordered array of geometric states (x, y, z, dx/dt, dy/dt, dz/dt,
    in kilometers and kilometers per second) of body relative to center, specified relative
    to frame. Useful for spice function "spkw09_c", for example.

    Parameters
    ----------
    ets : np.ndarray
        Array of TDB seconds from J2000 (et dates) of which the data will be taken
    pos_iau_earth : np.ndarray
        Rectangular coordinates of the point, referencing IAU_EARTH frame.
    delta_t : float
        TDB seconds between states
    source_frame : str
        Name of the frame to transform from.
    target_frame : str
        Name of the frame which the location will be referencing.

    Returns
    -------
    ndarray of float
        ndarray containing the states calculated
    """
    num_coordinates = 3
    n_state_attributes = 6
    states = np.zeros((len(ets), n_state_attributes))
    for i, et_value in enumerate(ets):
        states[i, :num_coordinates] = np.dot(
            spice.pxform(source_frame, target_frame, et_value),
            pos_iau_earth)

    for i in range(len(ets) - 1):
        states[i, num_coordinates:] = (states[i + 1, :num_coordinates] -
                                       states[i, :num_coordinates]) / delta_t

    pos_np1 = np.dot(
        spice.pxform(source_frame, target_frame, ets[-1] + delta_t),
        pos_iau_earth)
    states[-1, num_coordinates:] = (pos_np1 - states[-1, :num_coordinates]) / delta_t
    return states

@dataclass
class _EarthLocation():
    """
    Data for the creation of an observer point on earth surface

    Attributes
    ----------
    point_id : int
        ID code that will be associated with the point on Earth's surface
    states : np.ndarray of float64
        Array of geometric states of body relative to center
    """
    __slots__ = ['point_id', 'states']
    def __init__(self, point_id: int, lat: float, lon: float, altitude: float, ets: np.ndarray,
                 delta_t: float, source_frame: str, target_frame: str):
        """
        Parameters
        ----------
        point_id : int
            ID code that will be associated with the point on Earth's surface
        lat : float
            Geographic latitude of the observer point
        lon : float
            Geographic longitude of the observer point
        altitude : float
            Altitude over the sea level in meters.
        ets : np.ndarray
            Array of TDB seconds from J2000 (et dates) of which the data will be taken
        delta_t : float
            TDB seconds between states
        source_frame : str
            Name of the frame to transform from.
        target_frame : str
            Name of the frame which the location will be referencing.
        """
        self.point_id = point_id
        eq_rad = 6378 # Earth equatorial radius
        pol_rad = 6357 # Earth polar radius
        alt_km = altitude/1000
        flattening = (eq_rad - pol_rad)/eq_rad
        pos_iau_earth = spice.pgrrec('EARTH', math.radians(lon), math.radians(lat), alt_km,
                                     eq_rad, flattening)
        self.states = _calculate_states(ets, pos_iau_earth, delta_t, source_frame, target_frame)

def _get_moon_data(utc_time: str) -> MoonData:
    """Calculation of the moon data for the given utc_time for the loaded observer

    Parameters
    ----------
    utc_time : str
        Time at which the ELI will be calculated, in a valid UTC DateTime format

    Returns
    -------
    MoonData
        Moon data obtained from SPICE toolbox
    """
    et_date = spice.str2et(utc_time)

    _, radios_luna = spice.bodvrd("MOON", "RADII", 3)
    m_eq_rad = radios_luna[0] # 1738.1 # Moon Equatorial Radius
    m_pol_rad = radios_luna[2] # 1736 # Moon polar radius
    flattening = (m_eq_rad-m_pol_rad)/m_eq_rad

    # Calculate moon phase angle
    spoint, _, _ = spice.subpnt("INTERCEPT/ELLIPSOID", "MOON", et_date, 'MOON_ME',
                                "NONE", "Observer")
    phase = spice.phaseq(et_date, "MOON", "SUN", "Observer", "NONE")
    phase = math.degrees(phase)

    # Calculate selenographic coordinates of the observer
    sel_lon, sel_lat, _ = spice.recpgr("MOON", spoint, m_eq_rad, flattening)
    sel_lon = math.degrees(sel_lon)
    sel_lat = math.degrees(sel_lat)

    # Calculate selenographic longitude of sun
    sun_spoint, _, _ = spice.subslr("INTERCEPT/ELLIPSOID", "MOON", et_date, 'MOON_ME',
                                    "NONE", "SUN")
    sel_lon_sun_rad, _, _ = spice.recpgr("MOON", sun_spoint, m_eq_rad, flattening)

    # Calculate the distance between observer and moon (KM)
    state, _ = spice.spkezr("MOON", et_date, "MOON_ME", "NONE", "Observer")
    distance_observer_moon = math.sqrt(state[0]**2 + state[1]**2 + state[2]**2)

    # Calculate the distance between sun and moon (AU)
    state, _ = spice.spkezr("MOON", et_date, "MOON_ME", "NONE", "SUN")
    distance_sun_moon = math.sqrt(state[0]**2 + state[1]**2 + state[2]**2)
    distance_sun_moon = spice.convrt(distance_sun_moon, "KM", "AU")

    limit_lat = 90
    if sel_lat > limit_lat:
        sel_lat -= limit_lat*2
    elif sel_lat < -limit_lat:
        sel_lat += limit_lat*2

    limit_lon = 180
    if sel_lon > limit_lon:
        sel_lon -= limit_lon*2
    elif sel_lon < -limit_lon:
        sel_lon += limit_lon*2

    moon_data = MoonData(distance_sun_moon, distance_observer_moon, sel_lon_sun_rad, sel_lat,
                         sel_lon,phase)
    return moon_data

def _get_moon_datas_id(utc_times: List[str], kernels_path: str,
                       observer_id: int) -> List[MoonData]:
    """Calculation of needed MoonDatas from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon.

    Parameters
    ----------
    utc_times : list of str
        Times at which the ELI will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    observer_id : int
        Observer's body ID

    Returns
    -------
    list of MoonData
        Moon data obtained from SPICE toolbox
    """
    kernels = ["moon_pa_de421_1900-2050.bpc", "moon_080317.tf",
               "pck00010.tpc", "naif0011.tls", "de421.bsp", "custom.bsp", "earth_assoc_itrf93.tf",
               "earth_latest_high_prec.bpc", "earth_070425_370426_predict.bpc"]

    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        spice.furnsh(k_path)

    spice.boddef("Observer", observer_id)
    moon_datas = []

    for utc_time in utc_times:
        moon_datas.append(_get_moon_data(utc_time))

    spice.kclear()

    return moon_datas

def _create_earth_point_kernel(utc_times: List[str], kernels_path: str, lat: int, lon: int,
                               altitude: float, id_code: int) -> None:
    """Creates a SPK custom kernel file containing the data of a point on Earth's surface

    Parameters
    ----------
    utc_times : list of str
        Times at which the ELI will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    lat : float
        Geographic latitude (in degrees) of the location.
    lon : float
        Geographic longitude (in degrees) of the location.
    altitude : float
        Altitude over the sea level in meters.
    id_code : int
        ID code that will be associated with the point on Earth's surface
    """
    kernels = ["pck00010.tpc", "naif0011.tls", "earth_assoc_itrf93.tf",
               "de421.bsp", "earth_latest_high_prec.bpc", "earth_070425_370426_predict.bpc"]
    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        spice.furnsh(k_path)

    polynomial_degree = 5
    # Degree of the lagrange polynomials that will be used to interpolate the states
    delta_t = 1 # TDB seconds between states. Arbitrary.
    min_states_polynomial = polynomial_degree + 1
    # Min # states that are required to define a polynomial of that degree
    ets = np.array([])
    left_states = int(min_states_polynomial/2)
    right_states = left_states + min_states_polynomial%2
    for utc_time in utc_times:
        et0 = spice.str2et(utc_time)
        etprev = et0 - delta_t * left_states
        etf = et0 + delta_t * right_states
        ets_t = np.arange(etprev, etf, delta_t)
        for et_t in ets_t:
            if et_t not in ets:
                index = np.searchsorted(ets, et_t)
                ets = np.insert(ets, index, et_t)

    target_frame = source_frame = 'ITRF93'
    obs = _EarthLocation(id_code, lat, lon, altitude, ets, delta_t, source_frame, target_frame)

    custom_kernel_path = os.path.join(kernels_path, CUSTOM_KERNEL_NAME)
    handle = spice.spkopn(custom_kernel_path, 'SPK_file', 0)

    center = EARTH_ID_CODE
    spice.spkw09(handle, obs.point_id, center, target_frame,
                 ets[0], ets[-1], '0', polynomial_degree, len(ets),
                 obs.states.tolist(), ets.tolist())
    spice.spkcls(handle)

    spice.kclear()

def _remove_custom_kernel_file(kernels_path: str) -> None:
    """Remove the custom SPK kernel file if it exists

    Parameters
    ----------
    kernels_path : str
        Path where the SPICE kernels are stored
    """
    custom_kernel_path = os.path.join(kernels_path, CUSTOM_KERNEL_NAME)
    if os.path.exists(custom_kernel_path):
        os.remove(custom_kernel_path)

def get_moon_datas(lat: float, lon: float, altitude: float, utc_times: List[str],
                   kernels_path: str) -> List[MoonData]:
    """Calculation of needed Moon data from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon.

    Parameters
    ----------
    lat : float
        Geographic latitude (in degrees) of the location.
    lon : float
        Geographic longitude (in degrees) of the location.
    altitude : float
        Altitude over the sea level in meters.
    utc_times : str
        Times at which the ELI will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    Returns
    -------
    list of MoonData
        Moon data obtained from SPICE toolbox
    """
    id_code = 399100
    _remove_custom_kernel_file(kernels_path)
    _create_earth_point_kernel(utc_times, kernels_path, lat, lon, altitude, id_code)
    return _get_moon_datas_id(utc_times, kernels_path, id_code)
