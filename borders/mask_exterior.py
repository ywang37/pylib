"""
Created on September 2, 2020

@author: Yi Wang
"""

import cartopy.io.shapereader as shpreader
from copy import deepcopy
from matplotlib.collections import PathCollection
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import os

from mylib.cartopy_plot import cartopy_plot, scatter

#
#------------------------------------------------------------------------------
#
def pixel_in_shape(pixel, shape):
    """ Find all pixels in the shape
    (Yi Wang, 09/02/2020)

    Parameters
    ----------
    pixel : 1-D or 2-D numpy array
        If it is a single pixel, pixel is 1-D array. 
        Example: np.array([lon, lat])
        IF is multiple pixels, pixels is (N, 2) demension array.
    shape : list
        Elements are (Longitude, Latitude)

    Returns
    -------
    out_dict : dict
        flag : bool element of 1-D array
        shape : list of (Longitude, Latitude)

    """

    Path = mpath.Path

    # interior shape
    interior_data = []
    interior_data.append( (Path.MOVETO, shape[0]) )
    for i in range(1, len(shape)-1):
        interior_data.append( (Path.LINETO, shape[i]) )
    interior_data.append( (Path.CLOSEPOLY, shape[-1]) )
    in_codes, in_verts = zip(*interior_data)
    interior = mpath.Path(in_verts, in_codes)

    dim = pixel.shape
    N   = len(dim)
    if  N == 1:
        flag = interior.contains_point( pixel )
    elif N == 2:
        flag = interior.contains_points( pixel )
    else:
        print('pixel_in_shape: ERROR')
        print('pixel dimension is {}'.format(N))
        print('It should be no larger than 2')

    out_dict = {}
    out_dict['flag'] = flag
    out_dict['shape'] = shape

    return out_dict
#
#------------------------------------------------------------------------------
#
def pixel_in_contiguous_US(pixel):
    """ Find all pixels in Contiguous US
    (Yi Wang, 09/02/2020)

    Parameters
    ----------
    pixel : 1-D or 2-D numpy array
        If it is a single pixel, pixel is 1-D array. 
        Example: np.array([lon, lat])
        IF is multiple pixels, pixels is (N, 2) demension array.

    Returns
    -------
    out_dict : dict
        flag : bool element of 1-D array
        shape : list of (Longitude, Latitude)

    """

    curr_dir = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + '/'
    filename = curr_dir + 'data/cb_2017_us_nation_5m/cb_2017_us_nation_5m'
    shp = shpreader.Reader(filename)
    geometries = shp.geometries()
    geometry = next(geometries)
    shape = list(zip(*geometry[84].exterior.coords.xy))

    out_dict = pixel_in_shape(pixel, shape)

    return out_dict
#
#------------------------------------------------------------------------------
#
def pixel_in_US_state(pixel, state):
    """ Find all pixels in the specific state.
    (Yi Wang, 09/02/2020)

    Parameters
    ----------
    pixel : 1-D or 2-D numpy array
        If it is a single pixel, pixel is 1-D array. 
        Example: np.array([lon, lat])
        IF is multiple pixels, pixels is (N, 2) demension array.
    state : str
        US state

    Returns
    -------
    out_dict : dict
        flag : bool element of 1-D array
        shape : list of (Longitude, Latitude)

    """

    curr_dir = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + '/'
    filename = curr_dir + 'data/cb_2017_us_state_5m/cb_2017_us_state_5m'
    shp = shpreader.Reader(filename)

    # get geometry
    records = shp.records()
    state_list = []
    for rec in records:
        state_list.append(rec.attributes['NAME'])
        if (state == rec.attributes['NAME']):
            geometry = rec.geometry
    state_list.sort()
    if state not in state_list:
        print('state is ' + state + ', but it is incorrect.')
        print('Only state names blelow can be used.')
        for name in state_list:
            print(name)
        exit()

    # A geometry is a MultiPolygon. A MultiPolygon consits of
    # one or more than one Polygon; we only need the largest
    # Polygon. The rest Polygons are islands.
    shape = list(zip(*geometry[0].exterior.coords.xy))
    if (len(geometry) > 1):
        for i in range(1, len(geometry)):
            tmp_shape = list(zip(*geometry[i].exterior.coords.xy))
            if (len(tmp_shape) > len(shape)):
                shape = tmp_shape

    out_dict = pixel_in_shape(pixel, shape)

    return out_dict
#
#------------------------------------------------------------------------------
#
# test example
if '__main__' == __name__:


    min_lat = 25.0
    max_lat = 50.0
    min_lon = -125.0
    max_lon = -65.0

    N_lat = 100
    N_lon = 240
    lat_e = np.linspace(min_lat, max_lat, N_lat+1)
    lat_c = (lat_e[0:-1] + lat_e[1:]) * 0.5
    lon_e = np.linspace(min_lon, max_lon, N_lon+1)
    lon_c = (lon_e[0:-1] + lon_e[1:]) * 0.5

    lon_c, lat_c = np.meshgrid(lon_c, lat_c)
    lon_e, lat_e = np.meshgrid(lon_e, lat_e)

    pixel = np.array([lon_c.flatten(), lat_c.flatten()]).T


    ###############################
    # test pixel_in_contiguous_US
    ###############################

    # multiple pixels
    out_dict = pixel_in_contiguous_US(pixel)
    flag = np.reshape(out_dict['flag'], lat_c.shape)
    pout = cartopy_plot(lon_e, lat_e, flag, cbar=False,
            region_limit=[min_lat, min_lon, max_lat, max_lon],
            xtick=np.arange(min_lon, max_lon + 0.1, 20.0),
            ytick=np.arange(min_lat, max_lat + 0.1, 10.0),
            cl_res='10m', countries=True
            )

    # One pixel
    lat_in,  lon_in  = 40.0, -100.0
    lat_out, lon_out = 30.0, -70.0

    out_dict  = pixel_in_contiguous_US(np.array([lon_in,  lat_in ]))
    flag_in = out_dict['flag']
    out_dict = pixel_in_contiguous_US(np.array([lon_out, lat_out]))
    flag_out = out_dict['flag']
    print('Is lat = {}, lon = {} in the US? {}'.format(lat_in, 
        lon_in, flag_in))
    print('Is lat = {}, lon = {} in the US? {}'.format(lat_out, 
        lon_out, flag_out))
    ax = pout['ax']
    scatter(ax, lon_in, lat_in, color='blue')
    scatter(ax, lon_out, lat_out, color='lime')

    plt.savefig('pixel_in_contiguous_US.png', format='png', dpi=300)

    ###############################
    # pixel_in_US_state
    ##############################

    # multiple pixels
    out_dict = pixel_in_US_state(pixel, 'California')
    flag = np.reshape(out_dict['flag'], lat_c.shape)
    pout = cartopy_plot(lon_e, lat_e, flag, cbar=False,
            region_limit=[min_lat, min_lon, max_lat, max_lon],
            xtick=np.arange(min_lon, max_lon + 0.1, 20.0),
            ytick=np.arange(min_lat, max_lat + 0.1, 10.0),
            cl_res='10m', countries=True
            )

    plt.savefig('pixel_in_US_state.png', format='png', dpi=300)
#
#------------------------------------------------------------------------------
#
