"""
Created on September 17, 2019

@author: Yi Wang
"""

import numpy as np

def generate_grid(nlat_c, nlon_c):
    """ Generate latitudes and longitudes of grid edges and centers

    Parameters
    ----------
    n_lat_c : int
        # of grid box in the sourh-north direction
    n_lon_c : int
        # of grid box in the west-east direction

    Returns
    -------
    lat_e : 1D numpy array
        latidue edges
    lon_e : 1D numpy array
        longitude edges
    lat_c : 1D numpy array
        latitude edges
    lon_c : 1D numpy array
        longitude edges

    """

    lat_e = np.linspace( -90.0,  90.0, num=nlat_c+1)
    lon_e = np.linspace(-180.0, 180.0, num=nlon_c+1)

    lat_c = (lat_e[0:-1] + lat_e[1:]) * 0.5
    lon_c = (lon_e[0:-1] + lon_e[1:]) * 0.5

    return lat_e, lon_e, lat_c, lon_c

def get_center_index(edges, value):
    """ Get index of specific latitude or longitude.

    Parameters
    ----------
    edges : 1D numpy array
        latitude or longitude edges
    value : float
        latitude or longitude

    Returns
    -------
    i : int
       index
    
    """
    if (value < edges[0]) or (value > edges[-1]):
        print(' - get_center_index: !!! ERROR !!!, \
latitude or longitude is out of range')
        print(' latitide or longitude is {}'.format(value))
        print(' latitude or longitude edges are ', edges)
        exit()

    for i in range(len(edges)-1):
        if (value >= edges[i]) and (value <= edges[i+1]):
            return i

def get_center_index_2(value, start, end, step):
    """ Get index of specific latitude or longitude.

    Parameters
    ----------
    value : float
        latitude or longitude
    start : float
        Start of edge
    end : float
        End of edge
    step : float
        grid interval
    """

    if (value < start) or (value > end):
        print('value = {}'.format(value))
        print('start = {}'.format(start))
        print('end = {}'.format(end))
        print(' - get_center_index_2: out of range')

    N = round( (end - start) / step)

    i = int( (value - start) / step )
    if (i == N):
        i = i - 1

    return i

def region_limit_flag(lat, lon, region_limit):
    """ Find pixels in the *region_limit_flag*
    (Yi Wang, 11/27/2019)

    Parameters
    ----------
    lat : numpy array
        Latitude array
    lon : numpy array
        Longitude array
    region_limit : list-like
        [lat_min, lon_min, lat_max, lon_max]

    Returns
    -------
    flag : numpy logical array
        Mark pixels in the *region_limit_flag*

    """

    lat_min = region_limit[0]
    lon_min = region_limit[1]
    lat_max = region_limit[2]
    lon_max = region_limit[3]

    flag_lat = np.logical_and(lat >= lat_min, lat <= lat_max)
    flag_lon = np.logical_and(lon >= lon_min, lon <= lon_max)
    flag = np.logical_and(flag_lat, flag_lon)

    return flag

def generate_grid_gc_2x25():
    """
    """

    # resolution
    dlat = 2.0
    dlon = 2.5

    # latitude edge
    lat_e     = np.arange(  -91.0,   91.0001, dlat)
    lat_e[0]  = -90.0
    lat_e[-1] = 90.0

    # longitude edge
    lon_e = np.arange(-181.25, 178.7501, dlon)

    # latitude center
    lat_c = (lat_e[0:-1] + lat_e[1:]) / 2.0

    # longitude center
    lon_c = (lon_e[0:-1] + lon_e[1:]) / 2.0

    return lat_e, lon_e, lat_c, lon_c

def get_index_gc_2x25(lat, lon, lat_e, lon_e):
    """
    """

    # latitude
    if (lat < -90.0) or (lat > 90.0):
        print('get_index_2x25 error: latitue = {} is out of range.'.format(lat))
        exit()

    i = np.sum(lat_e <= lat) - 1
    if lat == 90.0:
        i = 90

    # longitude
    if (lon < -180.0) or (lon > 180.0):
        print('get_index_2x25 error: longitude = {} is out of range.'.format(lon))
        exit()

    if (lon <= 180.0) and (lon >= 178.75):
        j = 0
    else:
        j = np.sum(lon_e <= lon) - 1

    return i, j
