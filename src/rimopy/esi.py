"""ESI Extraterrestrial Solar Irradiation

This module contains the functionalities for obtaining the extraterrestrial solar irradiation
of a concrete wavelength, based on Wehrli (1985).

It exports the foollowing functions:

    * getESI - returns the expected extraterrestrial solar irradiation of a concrete wavelength
"""

import csv
from io import StringIO
from typing import Any, Tuple, Dict
import pkgutil
from scipy.interpolate import interp1d

def _getWehrliData() -> Dict[float, Tuple[float, float]]:
	"""Returns all wehrli data

	Returns
	-------
	A dict that has the wavelengths as keys (float), and as values it has tuples of the (W/sm/nm, W/sm) values.
	"""
	wehrli_bytes = pkgutil.get_data(__name__, 'data/wehrli.csv')
	wehrli_string = wehrli_bytes.decode()
	file = StringIO(wehrli_string)
	csvreader = csv.reader(file)
	next(csvreader) # Discard the header
	data = {}
	for row in csvreader:
		data[float(row[0])] = (float(row[1]), float(row[2]))
	file.close()
	return data

def _getClosestKey(d: Dict[float, Any], search_key: float) -> float:
	"""Finds the closest key of the dictionary

	Parameters
	----------
	d : dict of float and anything
		The dictionary from which we will find the closest key
	search_key : float
		The key of which we will find the closest one in the dictionary

	Returns
	-------
	float
		The closest key in the dictionary 
	"""
	return min(d.keys(), key = lambda key: abs(key-search_key))

def getESI(wavelength_nm: float) -> float:
	"""Gets the expected extraterrestrial solar irradiance at a concrete wavelength
    
    Parameters
    ----------
    wavelength_nm : float
        Wavelength (in nanometers) of which the extraterrestrial solar irradiance will be obtained
    
    Returns
    -------
    float
        The expected extraterrestrial solar irradiance in W/sm
	"""
	wehrli_data = _getWehrliData()
	wehrli_x = list(wehrli_data.keys())
	if wavelength_nm in wehrli_x:
		return wehrli_data[wavelength_nm][1]
	if wavelength_nm < wehrli_x[0]:
		return wehrli_data[0][1]
	if wavelength_nm > wehrli_x[-1]:
		return wehrli_data[-1][1]
	wehrli_y = list(map(lambda x : x[1], wehrli_data.values()))
	f = interp1d(wehrli_x, wehrli_y, 'cubic') # This works because, supposedly, python dicts preserve insertion order since 3.7
	return f(wavelength_nm).item()