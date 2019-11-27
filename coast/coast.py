"""
Created on November 26, 2019

@author: Yi Wang
"""

#from netCDF4 import Dataset
#import numpy as np

#
#------------------------------------------------------------------------------
#
def get_dist2coast(lat, lon, dist2coast, lat_st=-90.0, lon_st=-180.0,
        lat_int=0.04, lon_int=0.04):
    """ Given *lat* and *lon*, get distance [km]
    of the point to the closest coastline

    Parameters
    ----------
    lat : float
        Latitude of point
    lon : float
        Longitude of point
    dist2coast : 2-D numpy array 
        Distance to coastline (units: km)
    lat_st : float
        Starting latitude edge
    lon_st : float
        Starting longitude edge
    lat_int : float
        Latitude interval
    lon_int : float
        Longitude interval

    Returns
    -------
    dist : float
        Distance to the closest coastline [km]

    """

    # latitude and longitude indexes
    i = int( (lat - lat_st) / lat_int )
    j = int( (lon - lon_st) / lon_int )

    if (i == dist2coast.shape[0]):
        i = i - 1

    if (j == dist2coast.shape[1]):
        j = j -1

    if ( (i < 0) or (i >= dist2coast.shape[0]) ):
        print('i = {}'.format(i))
        print('lat = {}'.format(lat))
        print(' - get_dist2coast: latitude is out of range.')

    if ( (j < 0) or (j >= dist2coast.shape[1]) ):
        print('j = {}'.format(j))
        print('lon = {}'.format(lon))
        print(' - get_dist2coast: longitude is out of range.')

    dist = dist2coast[i,j]

    return dist
#
#------------------------------------------------------------------------------
#
