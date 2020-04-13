"""
Created on April 12, 2020

@author: Yi Wang
"""

import cartopy.io.shapereader as shpreader
import numpy as np
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

#
#------------------------------------------------------------------------------
#
def is_land(lon, lat):
    """ Determine is pixels are land.
    (ywang, 04/12/20)

    Parameters
    ----------
    lon : 2D numpy array
        Longitude
    lat : 2D numpy array
        Latitude

    Returns
    -------
    land_flag : 2D numpy array (bool)
        True means land

    """

    land_shp_fname = shpreader.natural_earth(resolution='10m',
            category='physical', name='land')

    land_geom = \
            unary_union(list(shpreader.Reader(land_shp_fname).geometries()))
    land = prep(land_geom)

    # function for one pixel
    def is_land_pixel(x, y):
        return land.contains(sgeom.Point(x, y))

    land_flag = np.full_like(lon, False)
    for i in range(lon.shape[0]):
        for j in range(lon.shape[1]):
            flag = is_land_pixel(lon[i,j], lat[i,j])
            land_flag[i,j] = flag

    return land_flag
#
#------------------------------------------------------------------------------
#
