"""Wehrli Gauss

This module converts the data from a csv file containing Wehrli (1985) data, to another csv
now containing the same data but through interpolation and a Gaussian filter.

This module exports the following functions:
    * filter_file - Converts de data from a file to another one.
"""

from scipy import ndimage
from scipy.interpolate import interp1d
import numpy as np
import csv
from typing import List

_INTERPOLATION_TYPE = 'linear'

def _interpolate_list(f: interp1d, minimum: float, maximum: float, step: float) -> List[float]:
    """Creates a list from an interpolation function, which elements will be from minimum
    to maximum, increasing step by step.

    Parameters:
    -----------
    f : interp1d
        Interpolation function
    minimum : float
        First element from which it will be generated the first interpolated value of the list
    maximum : float
        First element from which it might be generated the last interpolated value of the list.
        Depending on the value of step, it might not be reached.
    step: float
        The values generating the interpolated values will go from minimum until maximum, increaseng step by step.

    Returns:
    --------
    list of float
        List of the new interpolated values.
    """
    l = []
    for i in np.arange(minimum, maximum, step):
        l.append(f(i).item())
    return l

def _gaussian_filter_list(l: List[float], sigma: np.ScalarType) -> List[float]:
    """Gaussian filters the list

    Parameters:
    -----------
    l : list of float
        List of float values that will be filtered.
    sigma : scalar
        standard deviation for Gaussian kernel

    Returns:
    --------
    list of float
        Filtered list of values.
    """
    return ndimage.gaussian_filter1d(np.float_(l), sigma)

def _generate_new_list(data: List[List[float]]) -> List[List[float]]:
    """Generates the new interpolated and filtered list from the original list.

    Parameters:
    -----------
    data : list of list of float
        Original Wehrli data, obtained from the csv. Each element is another list, [nm, W/sm/nm, W/sm].

    Returns:
    --------
    list of list of float
        Interpolated and filtered list, following the same structure as the 'data' parameter.
    """
    f1 = interp1d([elem[0] for elem in data], [elem[1] for elem in data], _INTERPOLATION_TYPE)
    f2 = interp1d([elem[0] for elem in data], [elem[2] for elem in data], _INTERPOLATION_TYPE)
    minimum = 199.5
    maximum = 3000 #10075.0
    step = 1
    l1 = _interpolate_list(f1, minimum, maximum, step)
    l2 = _interpolate_list(f2, minimum, maximum, step)
    sigma = 1
    l1 = _gaussian_filter_list(l1, sigma)
    l2 = _gaussian_filter_list(l2, sigma)

    new_data = []
    index = 0
    for i in np.arange(minimum, maximum, step):
        new_elem = [i, l1[index], l2[index]]
        new_data.append(new_elem)
        index = index + 1
    return new_data

def filter_file(input: str, output: str):
    """Reads the data from 'input' path, interpolates and filters it, and then writes it in a file to 'output' path.
    
    Parameters:
    -----------
    input : str
        File path of the input data file.
    output : str
        Desired file path of the output file.
    """
    file = open(input, "r")
    csvreader = csv.reader(file)
    header = next(csvreader)
    data = []
    for row in csvreader:
        data.append([float(row[0]), float(row[1]), float(row[2])])
    file.close()
    data = _generate_new_list(data)

    file = open(output, "w")
    csvwriter = csv.writer(file)
    csvwriter.writerow(header)
    for elem in data:
        csvwriter.writerow(elem)    
    return data