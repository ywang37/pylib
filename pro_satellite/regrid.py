"""
Created on April 10, 2020

@author: Yi Wang
"""

import numpy as np

from mylib.grid_utility import generate_grid_2, get_2D_center_index_2

#
#------------------------------------------------------------------------------
#
def drop_in_the_box(grid_dict, sat_lat, sat_lon, sat_var,
        sat_flag=None):
    """ Regrid satellite data through drop in the box.
    (ywang, 04/10/20)

    Parameters
    ----------
    grid_dict : dict
        'proj' is a compulsory key. Its vaule can be
        'PlateCarree', ......
        'PlateCarree': Starting point, Ending point, and Step
            are given. Keys are 'lat_start', 'lat_end',
            lat_step', 'lon_start', 'lon_end' and 'lon_step'.
    sat_lat : 2-D numpy array
        Satellite latitude
    sat_lon : 2-D numpy array
        Satellite longitude
    sat_var : numpy array
        Satellite variable array to be (re)grided. It has
        at least 2 dimensions like (latitude, longitude, ...)
    sat_flag : logical 2-D numpy array or None (default)
        Determine whether a satellite variable is availble
        for (re)griding or not. True: available. False: not.
        If it is None, all observations are available for
        (re)griding.

    Returns
    -------
    out_dict : dict
        keys
        'sat_grid_ave'
        'sat_grid_count'
        'Latitude'
        'Longitude'
        'Latitude_e'
        'Longitude_e'

    """

    # satellite data dimensions
    sat_n_lat = sat_lat.shape[0]
    sat_n_lon = sat_lat.shape[1]

    # sat_flag
    if sat_flag is None:
        sat_flag = np.full_like(sat_lat, True)

    # grid
    if (grid_dict['proj'].lower() == 'PlateCarree'.lower()):

        lat_start = grid_dict['lat_start']
        lat_end   = grid_dict['lat_end']
        lat_step  = grid_dict['lat_step']
        lon_start = grid_dict['lon_start']
        lon_end   = grid_dict['lon_end']
        lon_step  = grid_dict['lon_step']
        lat_e, lon_e, lat, lon = \
                generate_grid_2(lat_step, lat_start, lat_end,
                        lon_step, lon_start, lon_end)
        n_lat = len(lat)
        n_lon = len(lon)
        print('lat_e: ')
        print(lat_e)
        print('lon_e: ')
        print(lon_e)
        print('lat_c: ')
        print(lat)
        print('lon_c: ')
        print(lon)
        Longitude, Latitude = np.meshgrid(lon, lat)
        Longitude_e, Latitude_e = np.meshgrid(lon_e, lat_e)

        # map satellite pixels to grid box
        sat_grid_ii = np.zeros( (sat_n_lat, sat_n_lon), dtype=int )
        sat_grid_jj = np.zeros( (sat_n_lat, sat_n_lon), dtype=int )
        for sat_i in range(sat_n_lat):
            for sat_j in range(sat_n_lon):

                # skip invalid satellite data
                if (not sat_flag[sat_i,sat_j]):
                    continue

                # get index
                grid_i, grid_j = get_2D_center_index_2(
                        sat_lat[sat_i,sat_j], 
                        lat_start, lat_end, lat_step, 
                        sat_lon[sat_i,sat_j], 
                        lon_start, lon_end, lon_step)

                # save index
                sat_grid_ii[sat_i,sat_j] = grid_i
                sat_grid_jj[sat_i,sat_j] = grid_j


    else:
        print(' - drop_in_the_box: grid_dict[\'proj\'] ' + \
                'value is ' + grid_dict['proj'] + \
                ', which is in correct. It can be \'PlateCarree\', ...')

    # shape of array that is used to save satellite observations
    # to new grid
    sat_grid_shape = list(sat_var.shape)
    sat_grid_shape[0] = n_lat
    sat_grid_shape[1] = n_lon
    sat_grid_shape = tuple(sat_grid_shape)

    # array to save satellite data
    sat_grid = np.zeros(sat_grid_shape, dtype=float)

    # array to save satellite data count
    sat_grid_count = np.zeros((n_lat,n_lon), dtype=int)

    # loop to process satellite data
    for sat_i in range(sat_n_lat):
        for sat_j in range(sat_n_lon):

            # skip invalid satellite data
            if (not sat_flag[sat_i,sat_j]):
                continue

            grid_i = sat_grid_ii[sat_i,sat_j]
            grid_j = sat_grid_jj[sat_i,sat_j]

            # accumulate satellite data at new grid
            sat_grid[grid_i,grid_j] += sat_var[sat_i,sat_j]

            # number of satellite data in a grid
            sat_grid_count[grid_i,grid_j] += 1

    # average of satellite data at new grid

    # broadcast
    sat_grid_count_full = np.broadcast_to(sat_grid_count.T,
            sat_grid.shape[::-1])
    sat_grid_count_full = sat_grid_count_full.T

    # average
    sat_grid_flag = (sat_grid_count_full > 0)
    sat_grid[sat_grid_flag] /= sat_grid_count_full[sat_grid_flag]
    sat_grid[np.logical_not(sat_grid_flag)] = np.nan

    # output dictionary
    out_dict = {}
    out_dict['Latitude'] = Latitude
    out_dict['Longitude'] = Longitude
    out_dict['Latitude_e'] = Latitude_e
    out_dict['Longitude_e'] = Longitude_e
    out_dict['sat_grid'] = sat_grid
    out_dict['sat_grid_count'] = sat_grid_count

    return out_dict
#
#------------------------------------------------------------------------------
#
