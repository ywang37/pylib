"""
Created on January 08, 2020

@author: Yi Wang
"""

from netCDF4 import Dataset
import numpy as np
from tqdm import tqdm

from mylib.grid_utility import generate_grid_2, get_center_index_2
from mylib.interpolate.interpolate import linear_interp_weight

def sat_model_sample(mod_coord_dict, mod_TAI93, mod_var_dict,
        sat_lat, sat_lon, sat_TAI93, sat_obs_dict,
        mod_flag=None, sat_flag=None,
        ):
    """ Sample model results according satellite observations
    and regrid satellite observations to model grids.

    Parameters
    ----------
    mod_coord_dict : dict
        *coord_format* is a compulsory key. Its vaule can be
        'ses', ......
        'ses' means Starting point, Ending point, and Step
            are given. Keys are 'lat_start', 'lat_end',
            'lat_step', 'lon_start', 'lon_end' and 'lon_step'
    mod_TAI93 : 1-D numpy array
        Seconds since 1993-01-01Z00:00:00 for model.
    mod_var_dict : dict
        Every value is a model variable array to be sampled.
        It has at least 3 dimensions like
        (time, latitude, longitude, ...)
    sat_lat : 2-D numpy array
        Satellite latitude
    sat_lon : 2-D numpy array
        Satellite longitude
    sat_TAI93 : 2-D numpy array
        Seconds since 1993-01-01Z00:00:00 for satellite.
    sat_obs_dict :
        Every value is a satellite variable array to
        be sampled. It has at least 2 dimensions like
        (latitude, longitude, ...)
    mod_flag : logical 2-D numpy array or None (default)
        Determine whether a grid is available for sampling
        or not. True: available. False: not.
        If it is None, all grids are available for sampling.
    sat_flag : logical 2-D numpy array or None (default)
        Determine whether a satellite variable is availble
        for regriding or not. True: available. False: not.
        If it is None, all observations are available for
        regriding.

    Returns
    -------
    out_dict : dict
        keys
        'sat_grid_dict'
        'sat_grid_num_dict'
        'sat_1D_dict'
        'sat_grid_dict'
        'sat_grid_num_dict'
        'sat_1D_dict'

    """

    out_dict = {}

    # model variable name list
    mod_var_name_list = list(mod_var_dict.keys())
    mod_var_name_list.sort()

    # model grid
    if (mod_coord_dict['coord_format'] == 'ses'):
        mod_lat_start = mod_coord_dict['mod_lat_start']
        mod_lat_end   = mod_coord_dict['mod_lat_end']
        mod_lat_step  = mod_coord_dict['mod_lat_step']
        mod_lon_start = mod_coord_dict['mod_lon_start']
        mod_lon_end   = mod_coord_dict['mod_lon_end']
        mod_lon_step  = mod_coord_dict['mod_lon_step']
        mod_lat_e, mod_lon_e, mod_lat, mod_lon = \
                generate_grid_2(mod_lat_step, mod_lat_start, mod_lat_end,
                        mod_lon_step, mod_lon_start, mod_lon_end)
        mod_n_lat = len(mod_lat)
        mod_n_lon = len(mod_lon)
        print('lat_e: ')
        print(mod_lat_e)
        print('lon_e: ')
        print(mod_lon_e)
        print('lat_c: ')
        print(mod_lat)
        print('lon_c: ')
        print(mod_lon)

    else:
        print(' - sat_model_sample: mod_coord_dict[\'coord_format\'] ' + \
                'value is ' + mod_coord_dict['coord_format'] + \
                ', which is incorrect. It can be \'ses\', ...')
        exit()

    # check mod_TAI93
    for i in range(len(mod_TAI93)-1):
        if (mod_TAI93[i] >= mod_TAI93[i+1]):
            print(' - sat_model_sample: mod_TAI93 is not ascending.')
            print('mod_TAI93: ')
            print(mod_TAI93)
            print('mod_TAI93[{}]: {}'.format(i,mod_TAI93[i]))
            print('mod_TAI93[{}]: {}'.format(i+1,mod_TAI93[i+1]))
            exit()

    # check model dimensions
    for mod_var_name in mod_var_name_list:

        mod_var_shape = mod_var_dict[mod_var_name].shape

        if ( mod_var_shape[0] != mod_TAI93.shape[0] ):
            print(' - sat_model_sample: ' + mod_var_name + \
                    ' time is inconsistent ' + \
                    'with that in mod_TAI93.')
            print('mod_TAI93: ')
            print(mod_TAI93)
            print('mod_TAI93.shape[0]: {}'.format(mod_TAI93.shape[0]))
            print(mod_var_name + \
                    '.shape[0]: {}'.format(mod_var_shape[0]))
            exit()

        if ( (mod_var_shape[1] != mod_n_lat) or \
                (mod_var_shape[2] != mod_n_lon) ):
            print(' - sat_model_sample: ' + mod_var_name  + \
                    ' grid is inconsistent ' + \
                    'with the dimension of model variable.')
            print(mod_coord_dict)
            print('mod_n_lat: {}'.format(mod_n_lat))
            print('mod_n_lon: {}'.format(mod_n_lon))
            print('mod_var_shape: ')
            print(mod_var_shape)
            print('lat_e: ')
            print(mod_lat_e)
            print('lon_e: ')
            print(mod_lon_e)
            print('lat_c: ')
            print(mod_lat_c)
            print('lon_c:')
            print(mod_lon_c)
            exit()


    # satellite data dimensions
    sat_n_lat = sat_lat.shape[0]
    sat_n_lon = sat_lat.shape[1]

    # mod_flag
    if mod_flag is None:
        mod_flag = np.full((mod_n_lat, mod_n_lon), True)
        print(mod_flag.shape)

    # sat_flag
    if sat_flag is None:
        sat_flag = np.full_like(sat_lat, True)
        print(sat_flag.shape)

    # mod_TAI93
    mod_TAI93_start = mod_TAI93[0]
    mod_TAI93_end   = mod_TAI93[-1]
    mod_n_time      = len(mod_TAI93)
    mod_sec_step    = (mod_TAI93_end - mod_TAI93_start) / (mod_n_time - 1)
    print(mod_TAI93_start)
    print(mod_TAI93_end)
    print(mod_n_time)
    print(mod_sec_step)

    # arrays that are used to save resmampled model variables
    # in the model grid
    mod_grid_dict = {}
    mod_grid_num_dict = {}
    mod_1D_dict = {}
    for mod_var_name in mod_var_name_list:

        # shape
        mod_grid_shape = mod_var_dict[mod_var_name].shape[1:]

        # model data
        mod_grid = np.zeros(mod_grid_shape, dtype=float)
        mod_grid_dict[mod_var_name] = mod_grid

        # number of model data
        mod_grid_num = np.zeros(mod_grid_shape[0:2], dtype=int)
        mod_grid_num_dict[mod_var_name] = mod_grid_num

        # save model data like station data
        mod_1D_dict[mod_var_name] = []

    # satellite observation name list and arrays that are used
    # to save satellite observations to model grid
    sat_obs_name_list = list(sat_obs_dict.keys())
    sat_obs_name_list.sort()
    sat_grid_dict = {}
    sat_grid_num_dict = {}
    sat_1D_dict = {}
    for sat_obs_name in sat_obs_name_list:

        sat_obs = sat_obs_dict[sat_obs_name]

        # shape of array that is used to save satellite observations
        # to model grid
        sat_grid_shape = list(sat_obs.shape)
        sat_grid_shape[0] = mod_n_lat
        sat_grid_shape[1] = mod_n_lon
        sat_grid_shape = tuple(sat_grid_shape)

        # satellite data
        sat_grid = np.zeros(sat_grid_shape, dtype=float)
        sat_grid_dict[sat_obs_name] = sat_grid

        # number of satellite
        sat_grid_num = np.zeros(sat_grid_shape[0:2], dtype=int)
        sat_grid_num_dict[sat_obs_name] = sat_grid_num

        # save satellite data like station data 
        sat_1D_dict[sat_obs_name] = []

    # sample model results according satellite observations
    # and regrid satellite observations to model grids.
    lat_ind_1D = []
    lon_ind_1D = []
    for sat_i in tqdm(range(sat_n_lat)):
        for sat_j in range(sat_n_lon):
            
            # skip invalid satellite data
            if (not sat_flag[sat_i,sat_j]):
                continue

            # skip satellite data that are not in
            # the time range of model
            if ( (sat_TAI93[sat_i,sat_j] < mod_TAI93_start) or \
                    (sat_TAI93[sat_i,sat_j] > mod_TAI93_end  ) ):
                continue

            # skip satellite data that are not in
            # the model grid region
            if ( (sat_lat[sat_i,sat_j] < mod_lat_start) or \
                    (sat_lat[sat_i,sat_j] > mod_lat_end) or \
                    (sat_lon[sat_i,sat_j] < mod_lon_start) or \
                    (sat_lon[sat_i,sat_j] > mod_lon_end) ):
                continue

            # determine the indexes of mod_TAI93, so satelite time
            # is within the corresponding model time
            ind_1 = int( (sat_TAI93[sat_i,sat_j] - mod_TAI93_start) 
                    / mod_sec_step )
            if (ind_1 == mod_n_time):
                ind_1 = ind_1 - 1
            ind_2 = ind_1 + 1
            #print(ind_1, mod_TAI93[ind_1])
            #print(ind_2, mod_TAI93[ind_2])
            if not ((sat_TAI93[sat_i,sat_j] >= mod_TAI93[ind_1]) or
                    (sat_TAI93[sat_i,sat_j] <= mod_TAI93[ind_2])):
                print(' - sat_model_sample: error')
                print('mod_TAI93: ')
                print(mod_TAI93)
                print('sat_TAI93[{},{}]: {}'.format(sat_i,sat_j,
                    sat_TAI93[sat_i,sat_j]))
                print('mod_TAI93[{}]: {}'.format(ind_1, mod_TAI93[ind_1]))
                print('mod_TAI93[{}]: {}'.format(ind_2, mod_TAI93[ind_2]))
                exit()

            # indexed of satellite data in model grid
            mod_i = get_center_index_2(sat_lat[sat_i,sat_j], 
                    mod_lat_start, mod_lat_end, mod_lat_step)
            mod_j = get_center_index_2(sat_lon[sat_i,sat_j],
                    mod_lon_start, mod_lon_end, mod_lon_step)
            lat_ind_1D.append(mod_i)
            lon_ind_1D.append(mod_j)

            # map satellite data to model grid
            for sat_obs_name in sat_obs_name_list:

                sat_obs      = sat_obs_dict[sat_obs_name]
                sat_grid     = sat_grid_dict[sat_obs_name]
                sat_grid_num = sat_grid_num_dict[sat_obs_name]
                sat_1D       = sat_1D_dict[sat_obs_name]

                # accumulate satellite data at model grid
                sat_grid[mod_i,mod_j] += sat_obs[sat_i,sat_j]

                # count satellite data at model grid
                sat_grid_num[mod_i,mod_j] += 1

                # save satellite data like station data
                sat_1D.append(sat_obs[sat_i,sat_j])

            # resample model variables according to satellite
            # overpass time
            for mod_var_name in mod_var_name_list:

                mod_var      = mod_var_dict[mod_var_name]
                mod_grid     = mod_grid_dict[mod_var_name]
                mod_grid_num = mod_grid_num_dict[mod_var_name]
                mod_1D       = mod_1D_dict[mod_var_name]

                # The closest model variables that are before
                # and after satellite overpass time
                mod_var_1 = mod_var[ind_1,mod_i,mod_j]
                mod_var_2 = mod_var[ind_2,mod_i,mod_j]

                # linear interpolation weight
                x  = sat_TAI93[sat_i,sat_j]
                x1 = mod_TAI93[ind_1]
                x2 = mod_TAI93[ind_2]
                weight = linear_interp_weight(x, x1, x2)
                #print(x, x1, x2)
                #print(weight)

                # linear interpolation
                mod_var_interp = \
                        (mod_var_1 * weight[0]) + (mod_var_2 * weight[1])

                # accumulte resampled model variables at model grid
                mod_grid[mod_i,mod_j] += mod_var_interp

                # count resampled model variables at model grid
                mod_grid_num[mod_i,mod_j] += 1

                # save resampled model variables like station data
                mod_1D.append(mod_var_interp)

    # average of satellite data at model grid
    for sat_obs_name in sat_obs_name_list:

        sat_grid     = sat_grid_dict[sat_obs_name]
        sat_grid_num = sat_grid_num_dict[sat_obs_name]

        # broadcast
        sat_grid_num_full = np.broadcast_to(sat_grid_num.T,
                sat_grid.shape[::-1])
        sat_grid_num_full = sat_grid_num_full.T
        print('test')
        print()

        # average 
        sat_grid_flag = (sat_grid_num_full > 0)
        sat_grid[sat_grid_flag] /= sat_grid_num_full[sat_grid_flag]
        sat_grid[np.logical_not(sat_grid_flag)] = np.nan

        # just convert list to array
        sat_1D_dict[sat_obs_name] = np.array(sat_1D_dict[sat_obs_name])

    # average of resampled model variables at model grid
    for mod_var_name in mod_var_name_list:

        mod_grid     = mod_grid_dict[mod_var_name]
        mod_grid_num = mod_grid_num_dict[mod_var_name]

        # broadcast
        mod_grid_num_full = np.broadcast_to(mod_grid_num.T,
                mod_grid.shape[::-1])
        mod_grid_num_full = mod_grid_num_full.T

        # average 
        mod_grid_flag = (mod_grid_num_full > 0)
        mod_grid[mod_grid_flag] /= mod_grid_num_full[mod_grid_flag]
        mod_grid[np.logical_not(mod_grid_flag)] = np.nan

        # just convert list to array
        mod_1D_dict[mod_var_name] = np.array(mod_1D_dict[mod_var_name])


    # lat and lon index
    lat_ind_1D = np.array(lat_ind_1D)
    lon_ind_1D = np.array(lon_ind_1D)
    ind_1D_dict = {}
    ind_1D_dict['lat_ind_1D'] = lat_ind_1D
    ind_1D_dict['lon_ind_1D'] = lon_ind_1D


    # save data to dictionary
    out_dict['sat_grid_dict']     = sat_grid_dict
    out_dict['sat_grid_num_dict'] = sat_grid_num_dict
    out_dict['sat_1D_dict']       = sat_1D_dict
    out_dict['mod_grid_dict']     = mod_grid_dict
    out_dict['mod_grid_num_dict'] = mod_grid_num_dict
    out_dict['mod_1D_dict']       = mod_1D_dict
    out_dict['ind_1D_dict']       = ind_1D_dict

    return out_dict

def save_sat_model_sample(filename, data_dict, save_2D=True, save_1D=True,
        verbose=True):
    """ Save dict returned from sat_model_sample. Model latitude and
    longitude information is also added to the dict before call the
    function.

    Parameters
    ----------
    filename : str
        netCDF file to save data.
    data_dict : dict
        Dictionary returned from sat_model_sample.  Model latitude and
        longitude information is also added to the dictionary
    save_2D : logical (default True)
        Save data in grid
    save_1D : logical (default True)
        Save data as station data

    Returns
    -------
    No returns

    """

    # get data
    sat_grid_dict     = data_dict['sat_grid_dict']
    sat_grid_num_dict = data_dict['sat_grid_num_dict']
    sat_1D_dict       = data_dict['sat_1D_dict']
    mod_grid_dict     = data_dict['mod_grid_dict']
    mod_grid_num_dict = data_dict['mod_grid_num_dict']
    mod_1D_dict       = data_dict['mod_1D_dict']
    ind_1D_dict       = data_dict['ind_1D_dict']

    mod_var_name_list = list(mod_grid_dict.keys())
    mod_var_name_list.sort()
    sat_obs_name_list = list(sat_grid_dict.keys())
    sat_obs_name_list.sort()

    if verbose:
        print(' - save_sat_model_sample: output ' + filename)

    # open file
    nc_f = Dataset(filename, 'w')

    # find if model variables have vertical dimension
    n_mod_lev = 1
    for mod_var_name in mod_var_name_list:
        if mod_grid_dict[mod_var_name].ndim == 3:
            n_mod_lev = mod_grid_dict[mod_var_name].shape[2]
            break

    # find if satellite observation have vertical dimension
    n_sat_lev = 1
    for sat_obs_name in sat_obs_name_list:
        if sat_grid_dict[sat_obs_name].ndim == 3:
            n_sat_lev = sat_grid_dict[sat_obs_name].shape[2]

    # fine if satellite observations have vertical dimension

    # grid, _e means edge
    Latitude    = data_dict['Latitude']
    Longitude   = data_dict['Longitude']
    Latitude_e  = data_dict['Latitude_e']
    Longitude_e = data_dict['Longitude_e']

    # Dimensions of a netCDF file
    if save_2D:
        dim_lat = nc_f.createDimension('Latitude',  Latitude.shape[0])
        dim_lon = nc_f.createDimension('Longitude', Latitude.shape[1])
        dim_lat_e = nc_f.createDimension('Latitude_e',  Latitude_e.shape[0])
        dim_lon_e = nc_f.createDimension('Longitude_e', Latitude_e.shape[1])
    
    if save_1D:
        for sat_var_name in sat_1D_dict:
            n_grid = sat_1D_dict[sat_var_name].shape[0]
            break
        dim_1D = nc_f.createDimension('grid', n_grid)

    if n_mod_lev > 1:
        dim_mod_lev = nc_f.createDimension('mod_lev', n_mod_lev)

    if n_sat_lev > 1:
        dim_sat_lev = nc_f.createDimension('sat_lev', n_sat_lev)

    # create variables in a netCDF file
    if save_2D:
        # lat and lon
        Latitude_v = nc_f.createVariable('Latitude', 'f4', 
                ('Latitude', 'Longitude'))
        Longitude_v = nc_f.createVariable('Longitude', 'f4',
                ('Latitude', 'Longitude'))
        Latitude_e_v = nc_f.createVariable('Latitude_e', 'f4', 
                ('Latitude_e', 'Longitude_e'))
        Longitude_e_v = nc_f.createVariable('Longitude_e', 'f4', 
                ('Latitude_e', 'Longitude_e'))
        # model variables
        nc_var_mod_grid_dict     = {}
        nc_var_mod_grid_num_dict = {}
        for mod_var_name in mod_var_name_list:
            mod_grid     = mod_grid_dict[mod_var_name]
            mod_grid_num = mod_grid_num_dict[mod_var_name]
            if mod_grid.ndim == 2:
                nc_var_mod_grid = \
                        nc_f.createVariable('mod_'+mod_var_name, 'f4',
                                ('Latitude', 'Longitude'))
            elif mod_grid.ndim == 3:
                nc_var_mod_grid = \
                        nc_f.createVariable('mod_'+mod_var_name, 'f4',
                                ('Latitude', 'Longitude', 'mod_lev'))
            else:
                print(' - save_sat_model_sample: mod_grid variable ' + 
                        mod_var_name + 
                        'has {} dimensions.'.format(mod_grid.ndim))
                exit()
            nc_var_mod_grid_num = \
                    nc_f.createVariable('mod_'+mod_var_name+'_num', 'f4',
                            ('Latitude', 'Longitude'))
            nc_var_mod_grid_dict[mod_var_name]     = nc_var_mod_grid
            nc_var_mod_grid_num_dict[mod_var_name] = nc_var_mod_grid_num
        # satellite observations
        nc_var_sat_grid_dict     = {}
        nc_var_sat_grid_num_dict = {}
        for sat_obs_name in sat_obs_name_list:
            sat_grid     = sat_grid_dict[sat_obs_name]
            sat_grid_num = sat_grid_num_dict[sat_obs_name]
            if sat_grid.ndim == 2:
                nc_var_sat_grid = \
                        nc_f.createVariable('sat_'+sat_obs_name, 'f4',
                                ('Latitude', 'Longitude'))
            elif sat_grid.ndim == 3:
                nc_var_sat_grid = \
                        nc_f.createVariable('sat_'+sat_obs_name, 'f4',
                                ('Latitude', 'Longitude', 'sat_lev'))
            else:
                print(' - save_sat_model_sample: sat_grid variable ' + 
                        sat_obs_name + 
                        'has {} dimensions.'.format(sat_grid.ndim))
                exit()
            nc_var_sat_grid_num = \
                    nc_f.createVariable('sat_'+sat_obs_name+'_num', 'f4',
                            ('Latitude', 'Longitude'))
            nc_var_sat_grid_dict[sat_obs_name]     = nc_var_sat_grid
            nc_var_sat_grid_num_dict[sat_obs_name] = nc_var_sat_grid_num

    if save_1D:
        # lat and lon index (0-based)
        lat_ind_1D_v = nc_f.createVariable('lat_ind_1D', 'int', ('grid',))
        lon_ind_1D_v = nc_f.createVariable('lon_ind_1D', 'int', ('grid',))
        # model variables
        nc_var_mod_1D_dict = {}
        for mod_var_name in mod_var_name_list:
            mod_1D = mod_1D_dict[mod_var_name]
            if mod_1D.ndim == 1:
                nc_var_mod_1D = \
                        nc_f.createVariable('mod_1D_'+mod_var_name, 'f4',
                                ('grid',))
            elif mod_1D.ndim == 2:
                nc_var_mod_1D = \
                        nc_f.createVariable('mod_1D_'+mod_var_name, 'f4',
                                ('grid', 'mod_lev'))
            else:
                print(' - save_sat_model_sample: mod_1D variable ' +
                        mod_var_name +
                        'has {} dimensions.'.format(mod_1D.ndim))
                exit()
            nc_var_mod_1D_dict[mod_var_name] = nc_var_mod_1D
        # satellite observations
        nc_var_sat_1D_dict = {}
        for sat_obs_name in sat_obs_name_list:
            sat_1D = sat_1D_dict[sat_obs_name]
            if sat_1D.ndim == 1:
                nc_var_sat_1D = \
                        nc_f.createVariable('sat_1D_'+sat_obs_name, 'f4',
                                ('grid',))
            elif sat_1D.ndim == 2:
                nc_var_sat_1D = \
                        nc_f.createVariable('sat_1D_'+sat_obs_name, 'f4',
                                ('grid', 'sat_lev'))
            else:
                print(' - save_sat_model_sample: sat_1D variable ' +
                        sat_obs_name +
                        'has {} dimensions.'.format(sat_1D.ndim))
                exit()
            nc_var_sat_1D_dict[sat_obs_name] = nc_var_sat_1D

    # write variables
    if save_2D:
        # lat and lon
        Latitude_v[:]    = Latitude
        Longitude_v[:]   = Longitude
        Latitude_e_v[:]  = Latitude_e
        Longitude_e_v[:] = Longitude_e
        # model variables
        for mod_var_name in mod_var_name_list:
            # value
            nc_var_mod_grid = nc_var_mod_grid_dict[mod_var_name]
            mod_grid = mod_grid_dict[mod_var_name]
            nc_var_mod_grid[:] = mod_grid
            # number
            nc_var_mod_grid_num = nc_var_mod_grid_num_dict[mod_var_name]
            mod_grid_num = mod_grid_num_dict[mod_var_name]
            nc_var_mod_grid_num[:] = mod_grid_num
        # satellite observations
        for sat_obs_name in sat_obs_name_list:
            # value
            nc_var_sat_grid = nc_var_sat_grid_dict[sat_obs_name]
            sat_grid = sat_grid_dict[sat_obs_name]
            nc_var_sat_grid[:] = sat_grid
            # number
            nc_var_sat_grid_num = nc_var_sat_grid_num_dict[sat_obs_name]
            sat_grid_num = sat_grid_num_dict[sat_obs_name]
            nc_var_sat_grid_num[:] = sat_grid_num

    if save_1D:
        # lat and lon
        lat_ind_1D_v = ind_1D_dict['lat_ind_1D']
        lon_ind_1D_v = ind_1D_dict['lon_ind_1D']
        # model variables
        for mod_var_name in mod_var_name_list:
            nc_var_mod_1D = nc_var_mod_1D_dict[mod_var_name]
            mod_1D = mod_1D_dict[mod_var_name]
            nc_var_mod_1D[:] = mod_1D
        # satellite observations
        for sat_obs_name in sat_obs_name_list:
            nc_var_sat_1D = nc_var_sat_1D_dict[sat_obs_name]
            sat_1D = sat_1D_dict[sat_obs_name]
            nc_var_sat_1D[:] = sat_1D

    # close file
    nc_f.close()
