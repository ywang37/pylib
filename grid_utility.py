"""
Created on September 17, 2019

@author: Yi Wang
"""

from area import area
from copy import deepcopy
import matplotlib
import numpy as np

from mylib.cartopy_plot import cartopy_plot

#
#------------------------------------------------------------------------------
#
def generate_grid(nlat_c, nlon_c, 
        lat_min=-90.0, lat_max=90.0, lon_min=-180.0, lon_max=180.0):
    """ Generate latitudes and longitudes of grid edges and centers

    Parameters
    ----------
    n_lat_c : int
        # of grid box in the sourh-north direction
    n_lon_c : int
        # of grid box in the west-east direction
    lat_min : float
        minimum latitude
    lat_max : float
        maximum latitude
    lon_min : float
        minimum longitude
    lon_max : float
        maximum longitude

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

    lat_e = np.linspace(lat_min, lat_max, num=nlat_c+1)
    lon_e = np.linspace(lon_min, lon_max, num=nlon_c+1)

    lat_c = (lat_e[0:-1] + lat_e[1:]) * 0.5
    lon_c = (lon_e[0:-1] + lon_e[1:]) * 0.5

    return lat_e, lon_e, lat_c, lon_c
#
#------------------------------------------------------------------------------
#
def generate_grid_2(lat_res, lat_min, lat_max,
        lon_res, lon_min, lon_max):
    """ Generate latitudes and longitudes of grid edges and centers
    
    Parameters
    ----------
    lat_res : float
        latitude resolution
    lon_res : float
        longitude resolution
    lat_min : float
        minimum latitude
    lat_max : float
        maximum latitude
    lon_min : float
        minimum longitude
    lon_max : float
        maximum longitude

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

    nlat_c = round( (lat_max - lat_min) / lat_res)
    nlon_c = round( (lon_max - lon_min) / lon_res)

    lat_e, lon_e, lat_c, lon_c = \
            generate_grid(nlat_c, nlon_c,
                    lat_min=lat_min, lat_max=lat_max, 
                    lon_min=lon_min, lon_max=lon_max)

    return lat_e, lon_e, lat_c, lon_c
#
#------------------------------------------------------------------------------
#
def get_pressure_index(pres_edge, pressure, check_flag=True):
    """ Get index of specific pressure layer
    (Yi Wang, 01/29/2020)

    Parameters
    ----------
    pres_edge : 1D numpy array
        Pressure levels in descening order.
    pressure : float
        pressure
    check_flag : bool (defaul: True)
        Check if *pres_edge* is in descending order

    Returns
    -------
    i : int
        index

    """

    # check range
    if (pressure > pres_edge[0]) or (pressure < pres_edge[-1]):
        print(' - get_pressure_index: !!! ERROR !!!, ' + \
                'pressure is out of range.')
        print(' pressure is {}'.format(pressure))
        print(' pressure edges are ', pres_edge)
        exit()

    # check descending order
    for i in range(len(pres_edge)-1):
        if (pres_edge[i+1] >= pres_edge[i]):
            print(' - get_pressure_index: pres_edge is not in' + \
                    'descending order.')
            print('pres_edge is: ')
            print(pres_edge)
            exit()

    # get index
    for i in range(len(pres_edge)-1):
        if ( (pressure <= pres_edge[i]) and (pressure >= pres_edge[i+1]) ):
            return i
#
#------------------------------------------------------------------------------
#
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
        return None

    for i in range(len(edges)-1):
        if (value >= edges[i]) and (value <= edges[i+1]):
            return i
#
#------------------------------------------------------------------------------
#
def get_center_index_latlon(lat_e, lon_e, region_limit):
    """ Get indexes for *region_limit*
    (ywang, 03/31/20)

    Parameters
    ----------
    lat_e : 1D numpy array
        Latitude edges
    lon_e : 1D numpy array
        Longitude edges
    region_limit : list_like
        [lat_min, lon_min, lat_max, lon_max]

    Returns
    -------
    i1, i2, j1, j2 : all ints

    """

    # lat_min
    i1 = get_center_index(lat_e, region_limit[0])

    # lat_max
    i2 = get_center_index(lat_e, region_limit[2])

    # lon_min
    j1 = get_center_index(lon_e, region_limit[1])

    # lon_max
    j2 = get_center_index(lon_e, region_limit[3])

    return i1, i2, j1, j2
#
#------------------------------------------------------------------------------
#
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
        return None

    N = round( (end - start) / step)

    i = int( (value - start) / step )
    if (i == N):
        i = i - 1

    return i
#
#------------------------------------------------------------------------------
#
def get_2D_center_index(lat_value, lat_edge, lon_value, lon_edge):
    """ Get index of specific latitude or longitude through
    call get_center_index.
    (ywang, 06/27/20)

    """

    i = get_center_index(lat_edge, lat_value)
    j = get_center_index(lon_edge, lon_value)

    return i, j
#
#------------------------------------------------------------------------------
#
def get_2D_center_index_2(lat_value, lat_start, lat_end, lat_step,
        lon_value, lon_start, lon_end, lon_step):
    """ Get index of specific latitude or longitude through
    call get_center_index_2.
    (ywang, 04/10/20)
    """
    
    i = get_center_index_2(lat_value, lat_start, lat_end, lat_step)
    j = get_center_index_2(lon_value, lon_start, lon_end, lon_step)

    return i, j
#
#------------------------------------------------------------------------------
#
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
#
#------------------------------------------------------------------------------
#
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
#
#------------------------------------------------------------------------------
#
def generate_grid_gc_na_05x0625():
    """
    """

    lat_min =    9.75
    lat_max =   70.25
    lat_res =    0.5
    lon_min = -140.3125
    lon_max =  -39.6875
    lon_res =    0.625

    lat_e, lon_e, lat_c, lon_c = generate_grid_2(lat_res, lat_min, lat_max,
            lon_res, lon_min, lon_max)

    return lat_e, lon_e, lat_c, lon_c
#
#------------------------------------------------------------------------------
#
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
#
#------------------------------------------------------------------------------
#
def grid_area_1(lat_e, lon_int, nlon):
    """ Calculate grid areas
    (Yi Wang, 02/17/2020)

    Parameters
    ----------
    lat_e : 1-D array
        Latitude edges
    lon_int : float
        Longitude interval
    nlon : int
        Number of grids along longitude

    Returns
    -------
    area_2D : 2-D array
        Grid areas. Unit is m^2

    Notes
    -----
    (1) The function may not work over polar region.

    """

    # area along latitude
    area_1D = np.zeros((len(lat_e)-1,))
    for i in range(len(area_1D)):
        coordinates = [ [[0.0, lat_e[i]], 
                         [0.0, lat_e[i+1]], 
                         [lon_int, lat_e[i+1]],
                         [lon_int, lat_e[i]],
                         [0.0, lat_e[i]]
                         ] ]
        area_1D[i] = area( {
            'type': 'Polygon',
            'coordinates': coordinates
            } )

    area_2D = np.tile(area_1D, nlon)
    area_2D = np.reshape(area_2D, (nlon,len(lat_e)-1))
    area_2D = np.transpose(area_2D)

    return area_2D
#
#------------------------------------------------------------------------------
#
def region_ave_sum(in_data, in_weight=None,
        flag_mask=None, 
        lat=None, lon=None, region_limit=None):
    """ Calculate average and sum in a region.
    (Yi Wang, 06/30/2020)

    Parameters
    ----------
    in_data : 2-D array
        It should be a 2-D numpy
    in_weight : 2-D array or None
        If None, weights are 1 for every grid.
    flag_mask : 2-D bool array
        If None, all elements are False, which means no mask.
    lat : 2-D array or None
        Center latitudes
    lon : 2-D array or None
        Center longitudes
    region_limit : list-like or None
        [lat_min, lon_min, lat_max, lon_max]
        If lat, lon, and region_limit are not None,
        only use data in the region_limit.

    Returns
    -------
    out_dict : dictionary
        Keys: 'mean', 'sum', 'data', 'weight', 'final_flag'

    """

    # weight
    if in_weight is None:
        weight = np.full_like(in_data, 1.0, dtype=float)
    else:
        weight = deepcopy(in_weight)

    # data
    data = deepcopy(in_data)

    # data mask
    if isinstance(data, np.ma.core.MaskedArray):
        data_mask = data.mask
    else:
        data_mask = np.full(data.shape, False)

    # data nan
    flag_data_nan = np.isnan(data)

    # flag_mask
    if flag_mask is None:
        flag_mask = np.full(data.shape, False)

    # region_mask
    if (lat is not None) and (lat is not None) and \
            (region_limit is not None):
        region_mask= np.logical_not( \
                region_limit_flag(lat, lon, region_limit) )
    else:
        region_mask = np.full(data.shape, False)

    # final flag for elelemts that are not considered.
    final_flag = np.logical_or(region_mask, np.logical_or(data_mask, \
            np.logical_or(flag_data_nan, flag_mask)))

    # nan 
    data[final_flag]   = np.nan
    weight[final_flag] = np.nan

    # out_dict
    out_dict = {}

    # sum
    out_dict['sum'] = np.nansum(data * weight)

    # mean
    out_dict['mean'] = out_dict['sum'] / np.nansum(weight)

    # data
    out_dict['data'] = data

    # weight
    out_dict['weight'] = weight

    # final_flag
    out_dict['final_flag'] = final_flag

    return out_dict
#
#------------------------------------------------------------------------------
#
def fill_region_color(lon, lat, region_flag_list, color_list,
    region_limit=None, ax=None, fig=None,
    xtick=None, ytick=None, cl_res=None,
    countries=True, states=True):
    """ Fill regions by different colors.
    (ywang, 06/12/2020)

    Parameters
    ----------
    lon : 2-D numpy array
        2-D longitude
    lat : 2-D numpy array
        2-D latitude array
    region_flag_list : list
        Elements are region_flag, and a region_flag is a 2-D
        bool array with True means the girds need to be filled.
    color_list : list
        Filled colors for each region.
    region_limit : tuple-like or None
        (min_lat, min_lon, max_lat, max_lon)
    ax : GeoAxes or None (default)
        Create a GeoAxes if ax is None. 
    fig : plt.figure() or None (default)
        Create a plt.figure() if both fig and ax are None
    xtick : list-like
        Longitude ticks
    ytick : list-like
        Latitude ticks
    cl_res : str
        Coastline resolution. 
        Currently can be one of “110m”, “50m”, and “10m”

    Returns
    -------
    out_dict :
        keys : ax, fig, mesh

    """

    # region index
    region_ind = np.full_like(region_flag_list[0], np.nan, dtype=float)

    for i in range(len(region_flag_list)):

        region_flag = region_flag_list[i]

        region_ind[region_flag] = i + 0.5

    # plot
    out_dict = cartopy_plot(lon, lat, region_ind, ax=ax, fig=fig,
            region_limit=region_limit, xtick=xtick, ytick=ytick,
            cmap=matplotlib.colors.ListedColormap(color_list),
            bad_c='white',
            cbar=False, vmin=0.0, vmax=len(region_flag_list),
            cl_res=cl_res, countries=countries, states=states)

    return out_dict
#
#------------------------------------------------------------------------------
#
