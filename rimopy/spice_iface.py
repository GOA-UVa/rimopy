"""SPICE iface

Interface with NASA's SPICE toolbox

It exports the following functions:

    * get_moon_datas - Calculates needed MoonData from SPICE toolbox
    * get_moon_datas_from_extra_kernels - Calculates needed MoonData from SPICE toolbox
        and using data from extra kernels for the observer body
"""

from typing import List

import numpy as np
from spicedmoon import spicedmoon

from .types import MoonDatas


def get_moon_datas_from_extra_kernels(
    utc_times: List[str],
    kernels_path: str,
    extra_kernels: List[str],
    extra_kernels_path: str,
    observer_name: str,
) -> MoonDatas:
    """Calculation of needed Moon data from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon.

    Parameters
    ----------
    utc_times : str
        Times at which the ELI will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    extra_kernels: list of str
        Custom kernels from which the observer body will be loaded, instead of calculating it.
    extra_kernels_path: str
        Folder where the extra kernels are located.
    observer_name: str
        Name of the body of the observer that will be loaded from the extra kernels.
    Returns
    -------
    MoonDatas
        Moon data obtained from SPICE toolbox
    """
    spmds = spicedmoon.get_moon_datas_from_extra_kernels(
        utc_times,
        kernels_path,
        extra_kernels,
        extra_kernels_path,
        observer_name,
        "ITRF93",
        ignore_bodvrd=False,
    )
    mds = []
    for spmd in spmds:
        mds.append(
            [
                spmd.dist_sun_moon_au,
                spmd.dist_obs_moon,
                spmd.lon_sun_rad,
                spmd.lat_obs,
                spmd.lon_obs,
                spmd.mpa_deg,
            ]
        )
    mds = np.array(mds).T
    mds = MoonDatas(mds[0], mds[1], mds[2], mds[3], mds[4], mds[5])
    return mds


def get_moon_datas(
    lat: float, lon: float, altitude: float, utc_times: List[str], kernels_path: str
) -> MoonDatas:
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
    MoonDatas
        Moon data obtained from SPICE toolbox
    """
    spmds = spicedmoon.get_moon_datas(
        lat, lon, altitude, utc_times, kernels_path, ignore_bodvrd=False
    )
    mds = []
    for spmd in spmds:
        mds.append(
            [
                spmd.dist_sun_moon_au,
                spmd.dist_obs_moon,
                spmd.lon_sun_rad,
                spmd.lat_obs,
                spmd.lon_obs,
                spmd.mpa_deg,
            ]
        )
    mds = np.array(mds).T
    mds = MoonDatas(mds[0], mds[1], mds[2], mds[3], mds[4], mds[5])
    return mds
