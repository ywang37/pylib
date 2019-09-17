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


