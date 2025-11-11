"""
Geometry interface utilities

This module provides the internal helper function used to resolve or construct
`MoonDatas` objects from different geometry sources. It acts as a bridge between
the SPICE interface and the higher-level ELI/ELRef modules.

Exports
-------
Functions
    resolve_mds
        Unifies the different geometry inputs (precomputed MoonDatas, EarthPoint,
        or extra SPICE kernels) into a single `MoonDatas` instance for use in RIMO
        irradiance and reflectance calculations.
"""

from typing import Optional, Union, List

from .types import MoonDatas, EarthPoint
from . import spice_iface

def resolve_mds(
    mds: Optional[MoonDatas] = None,
    earth_data: Optional[EarthPoint] = None,
    kernels_path: Optional[str] = None,
    utc_times: Optional[Union[str, List[str]]] = None,
    extra_kernels: Optional[List[str]] = None,
    extra_kernels_path: Optional[str] = None,
    observer_name: Optional[str] = None,
) -> MoonDatas:
    """
    Resolve or construct a `MoonDatas` object from one of the supported geometry sources.

    This utility unifies the three possible ways to obtain lunar observation geometry for
    RIMO calculations:

      - **A)** Use a precomputed `MoonDatas` instance (`mds`).
      - **B)** Derive geometry from an Earth surface point (`earth_data`) using standard SPICE kernels.
      - **C)** Derive geometry from custom observer kernels (`extra_kernels`).

    Exactly one of these sources must be provided.  If multiple or none are specified, a
    `ValueError` is raised.

    Parameters
    ----------
    mds : MoonDatas, optional
        Precomputed Moon geometry and distances (source A).
    earth_data : EarthPoint, optional
        Geographic coordinates and times of the observation (source B).
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

    Returns
    -------
    MoonDatas
        Fully populated `MoonDatas` object with all geometry parameters required for
        ELI or reflectance calculations.

    Raises
    ------
    ValueError
        If none or more than one geometry source is specified.
    """
    source_flags = [
        mds is not None,
        earth_data is not None and kernels_path is not None,
        (utc_times is not None and kernels_path is not None
         and extra_kernels is not None and extra_kernels_path is not None
         and observer_name is not None),
    ]
    if sum(bool(f) for f in source_flags) != 1:
        raise ValueError(
            "Specify exactly one geometry source: "
            "(A) mds, or (B) earth+kernels_path, or "
            "(C) utc_times+kernels_path+extra_kernels+extra_kernels_path+observer_name."
        )

    if mds is not None:
        return mds
    if earth_data is not None:
        return spice_iface.get_moon_datas(
            earth_data.lat, earth_data.lon, earth_data.altitude, earth_data.utc_times, kernels_path
        )
    return spice_iface.get_moon_datas_from_extra_kernels(
        utc_times, kernels_path, extra_kernels, extra_kernels_path, observer_name
    )
