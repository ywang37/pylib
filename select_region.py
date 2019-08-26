"""
Created on August 26, 2019

@author: Yi Wang
"""

import numpy as np
from netCDF4 import Dataset

def is_in_the_region(lat, lon, region):
    """ Check if the pixels are in the region.

    Parameters
    ----------
    lat : array
        Latitude.
    lon : array
        Longitude.
    region : tuple-like
        Define the region.
        (min_lat, min_lon, max_lat, max)

    Returns
    -------
    flag : array
        The pixel is in the region (True) or not (Flase).

    """

    min_lat = region[0]
    min_lon = region[1]
    max_lat = region[2]
    max_lon = region[3]

    # Latitude limit
    flag_lat = np.logical_and(lat > min_lat, lat < max_lat)

    # Longtidue limit
    flag_lon = np.logical_and(lon > min_lon, lon < max_lon)

    # Latidue and longitude limits
    flag = np.logical_and(flag_lat, flag_lon)

    return flag

def MODIS_in_the_region(filename, region):
    """ Check if there is at least one pixel of the granule (filename)
    that is in the region.

    Parameters
    ----------
    filename : str
        MODIS filename
    region : tuple-like
        Define the region.
        (min_lat, min_lon, max_lat, max)

    Returns
    -------
    flag : logical

    """

    # Read latitude and longitude.
    f = Dataset(filename, 'r')
    lat = f.variables['Latitude'][:]
    lon = f.variables['Longitude'][:]
    f.close()

    # Check every pixel.
    all_pixel_flag = is_in_the_region(lat, lon, region)

    # At least one pixel.
    flag = np.any(all_pixel_flag)

    return flag

def all_files_in_the_region(all_files, region):
    """ Return all the filenames that has at least one pixel of the 
        granule (filename) that is in the region.

    Parameters
    ----------
    all_files : tuple-like
        Element is filename.
    region : tuple-like
        Define the region.
        (min_lat, min_lon, max_lat, max)

    Returns
    -------
    select_file : list
        All files that are in the region

    """

    select_file = []
    for filename in all_files:

        flag = MODIS_in_the_region(filename, region)

        if flag:
            select_file.append(filename)

    return select_file
