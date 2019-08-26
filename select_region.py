


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

