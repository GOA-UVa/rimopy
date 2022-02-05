"""SPICE iface

Interface with NASA's SPICE toolbox

It exports the following functions:
    
    * getMoonData - Calculates needed MoonData from SPICE toolbox
"""

from operator import eq
import spiceypy as spice
import os
from .MoonData import MoonData
import math
import numpy as np

CUSTOM_KERNEL_NAME = "custom.bsp"
EARTH_ID_CODE = 399

class _EarthLocation():
    """
    Data for the creation of an observer point on earth surface

    Attributes
    ----------
    point_id : int
        ID code that will be associated with the point on Earth's surface
    ets : np.ndarray
        Array of TDB seconds from J2000 (et dates) of which the data will be taken
    states : np.ndarray of float64
        Array of geometric states of body relative to center
    """
    __slots__ = ['point_id', 'ets', 'states']
    def __init__(self, point_id: int, lat: float, lon: float, altitude: float, ets: np.ndarray, delta_t: float, min_states_polynomial: int, frame: str):
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
        min_states_polynomial : int
            Minimum number states that are required to define a Lagrange polynomial of the degree it's going to be defined
        frame : str
            Name of the frame which the location will be referencing.
        """
        num_coordinates = 3
        n_state_attributes = 6
        self.point_id = point_id
        eq_rad = 6378 # Earth equatorial radius
        pol_rad = 6357 # Earth polar radius
        alt_km = altitude/1000
        flattening = (eq_rad - pol_rad)/eq_rad
        pos_iau_earth = spice.pgrrec( 'EARTH', math.radians(lon), math.radians(lat), alt_km, eq_rad, flattening)
        states = np.zeros( ( len( ets ), n_state_attributes ) )
        for n in range( len( ets ) ):
            states[ n, :num_coordinates ] = np.dot(
                spice.pxform( 'IAU_EARTH', frame, ets[ n ] ),
                pos_iau_earth )

        for n in range( len( ets ) - 1 ):
            states[ n, num_coordinates: ] = ( states[ n + 1, :num_coordinates ] - states[ n, :num_coordinates ] ) / delta_t

        pos_np1 = np.dot(
                spice.pxform( 'IAU_EARTH', frame, ets[ -1 ] + delta_t ),
                pos_iau_earth )
        states[ -1, num_coordinates: ] = ( pos_np1 - states[ -1, :num_coordinates ] ) / delta_t
        self.states = states

def _getMoonDataID(utc_time: str, kernels_path: str, id: int) -> MoonData:
    """Calculation of needed Moon data from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon. 

    Parameters
    ----------
    utc_time : str
        Time at which the ELI will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    id : int
        Observer's body ID

    Returns
    -------
    'MoonData'
        Moon data obtained from SPICE toolbox
    """
    kernels = ["moon_pa_de421_1900-2050.bpc", "moon_080317.tf", "moon_assoc_me.tf", "pck00010.tpc", "naif0012.tls", "de440.bsp", "custom.bsp"]
    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        spice.furnsh(k_path)

    spice.boddef("Observer", id)

    et_date = spice.str2et(utc_time)

    re = 1738.1 # Moon Equatorial Radius
    rp = 1736 # Moon polar radius
    f = (re-rp)/re

    # Calculate moon phase angle
    spoint, _, _ = spice.subpnt("NEAR POINT/ELLIPSOID", "MOON", et_date, 'MOON_ME', "NONE", "Observer")
    _, _, phase, _, _ = spice.ilumin("ELLIPSOID", "MOON", et_date, "MOON_ME", "NONE", "Observer", spoint)
    phase = math.degrees(phase)

    # Calculate selenographic coordinates of the observer
    sel_lon, sel_lat, _ = spice.recpgr("MOON", spoint, re, f)
    sel_lon = math.degrees(sel_lon)
    sel_lat = math.degrees(sel_lat)

    # Calculate selenographic longitude of sun
    sun_spoint, _, _ = spice.subpnt("NEAR POINT/ELLIPSOID", "MOON", et_date, 'IAU_MOON', "NONE", "SUN")
    sel_lon_sun_rad, _, _ = spice.recpgr("MOON", sun_spoint, re, f)

    # Calculate the distance between observer and moon (KM)
    obs_pos, _ = spice.spkpos("MOON", et_date, "MOON_ME", "NONE", "Observer")
    distance_observer_moon = spice.vnorm( obs_pos )

    # Calculate the distance between sun and moon (AU)
    sun_pos, _ = spice.spkpos("MOON", et_date, "J2000", "NONE", "SUN")
    distance_sun_moon = spice.vnorm( sun_pos )
    distance_sun_moon = spice.convrt( distance_sun_moon, "KM", "AU")    

    spice.kclear()

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

    md = MoonData(distance_sun_moon, distance_observer_moon, sel_lon_sun_rad, sel_lat, sel_lon, phase)

    return md

def _createEarthPointKernel(utc_time: str, kernels_path: str, lat: int, lon: int, altitude: float, id_code: int) -> None:
    """Creates a SPK custom kernel file containing the data of a point on Earth's surface

    Parameters
    ----------
    utc_time : str
        Time at which the ELI will be calculated, in a valid UTC DateTime format
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
    pck_kernel_path = os.path.join(kernels_path, "pck00010.tpc")
    naif_kernel_path = os.path.join(kernels_path, "naif0012.tls")
    spice.furnsh(pck_kernel_path)
    spice.furnsh(naif_kernel_path)

    polynomial_degree = 5 # Degree of the lagrange polynomials that will be used to interpolate the states
    delta_t = 1000 # TDB seconds between states
    et0 = spice.str2et(utc_time)
    min_states_polynomial = polynomial_degree + 1 # Min # states that are required to define a polynomial of that degree
    etf = et0 + delta_t * min_states_polynomial
    ets = np.arange(et0, etf, delta_t)

    frame = 'J2000'
    obs = _EarthLocation(id_code, lat, lon, altitude, ets, delta_t, min_states_polynomial, frame)

    custom_kernel_path = os.path.join(kernels_path, CUSTOM_KERNEL_NAME)
    handle = spice.spkopn(custom_kernel_path, 'SPK_file', 0 )
    
    center = EARTH_ID_CODE
    spice.spkw09( handle, obs.point_id, center, frame,
        ets[ 0 ], ets[ -1 ], '0', polynomial_degree, len( ets ),
        obs.states.tolist(), ets.tolist() )
    spice.spkcls(handle)

    spice.kclear()

def _removeCustomKernelFile(kernels_path: str) -> None:
    """Remove the custom SPK kernel file if it exists

    Parameters
    ----------
    kernels_path : str
        Path where the SPICE kernels are stored
    """
    custom_kernel_path = os.path.join(kernels_path, CUSTOM_KERNEL_NAME)
    if os.path.exists(custom_kernel_path):
        os.remove(custom_kernel_path)

def getMoonData(lat: float, long: float, altitude: float , utc_time: str, kernels_path: str) -> MoonData:
    """Calculation of needed Moon data from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon. 

    Parameters
    ----------
    lat : float
        Geographic latitude (in degrees) of the location.
    long : float
        Geographic longitude (in degrees) of the location.
    altitude : float
        Altitude over the sea level in meters.
    utc_time : str
        Time at which the ELI will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    Returns
    -------
    'MoonData'
        Moon data obtained from SPICE toolbox
    """
    id_code = 301100
    _removeCustomKernelFile(kernels_path)
    _createEarthPointKernel(utc_time, kernels_path, lat, long, altitude, id_code)
    return _getMoonDataID(utc_time, kernels_path, id_code)
    