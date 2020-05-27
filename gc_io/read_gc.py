"""
Created on May 26, 2020

@author: Yi Wang
"""

import datetime
from netCDF4 import Dataset
import numpy as np
import os

from mylib.io import read_nc

def read_inst_resample(file_list, varname_list, latlon_flag=True,
        time=True):
    """ read data from GC instant output files for
    resampling according to satellite overpass time.

    Parameters
    ----------
    file_list : list
        A list of instant output files. It usually includes
        files of previous one day, current day, and next
        one day.
    varname_list : list
        A list of variable names.
    latlon_flag : bool
        Whether or not get latitude and longitude
    time : bool
        Whether or not get time.

    Returns
    -------

    """

    out_dict = {}

    # 3-D or 4-D variables
    for varname in varname_list:
        out_dict[varname] = []

    if time:
        date_time = []

    # read data
    for i in range(len(file_list)):

        filename = file_list[i]

        # Determine if file exists
        if not os.path.exists(filename):
            print(' - read_inst_resample: WARNING! ' + filename + \
                    ' does not exist.')
            continue

        # get variables
        data_1 = read_nc(filename, varname_list, verbose=True)
        for varname in varname_list:
            out_dict[varname].append(data_1[varname])

        # get time
        if time:

            fid = Dataset(filename, 'r')

            time_var = fid.variables['time']

            time_start = getattr(time_var, 'units')
            time_start = time_start.split()
            dt_string = time_start[2] + ' ' + time_start[3]
            time_start = datetime.datetime.strptime(dt_string,
                    '%Y-%m-%d %H:%M:%S')

            time_delta = time_var[:]

            for i in range(len(time_delta)):
                date_time.append(time_start + \
                        datetime.timedelta(minutes=time_delta[i]))

            fid.close()


        # get latitude and longitude
        if latlon_flag:

            latlon_flag = False

            data_3 = read_nc(filename, ['lat', 'lon'], verbose=True)

            # centers
            out_dict['latitude']  = data_3['lat']
            out_dict['longitude'] = data_3['lon']

            # latitudue edges
            out_dict['latitude_e'] = \
                    ( data_3['lat'][:-1] + data_3['lat'][1:] ) * 0.5
            lat_del_1 = out_dict['latitude_e'][1] - out_dict['latitude_e'][0]
            lat_del_2 = out_dict['latitude_e'][2] - out_dict['latitude_e'][1]
            if ((lat_del_2 - lat_del_1) < 1e-6):
                out_dict['latitude_e'] = np.insert(out_dict['latitude_e'],
                        0, out_dict['latitude_e'][0] - lat_del_2)
                out_dict['latitude_e'] = np.insert(out_dict['latitude_e'],
                        len(out_dict['latitude_e']),
                        out_dict['latitude_e'][-1] + lat_del_2)
            else:
                out_dict['latitude_e'][0] = \
                        out_dict['latitude_e'][1] - lat_del_2
                out_dict['latitude_e'][-1] = \
                        out_dict['latitude_e'][-2] + lat_del_2
                lat_del_half = lat_del_2 * 0.5
                out_dict['latitude_e'] = np.insert(out_dict['latitude_e'],
                        0, out_dict['latitude_e'][0] - lat_del_half)
                out_dict['latitude_e'] = np.insert(out_dict['latitude_e'],
                        len(out_dict['latitude_e']),
                        out_dict['latitude_e'][-1] + lat_del_half)


            # longitude edges
            out_dict['longitude_e'] = \
                    ( data_3['lon'][:-1] + data_3['lon'][1:] ) * 0.5
            lon_del = out_dict['longitude_e'][1] - out_dict['longitude_e'][0]
            out_dict['longitude_e'] = np.insert(out_dict['longitude_e'],
                    0, out_dict['longitude_e'][0] - lon_del)
            out_dict['longitude_e'] = np.insert(out_dict['longitude_e'],
                    len(out_dict['longitude_e']),
                    out_dict['longitude_e'][-1] + lon_del)

    # conver date_time to TAI93
    if time:

        out_dict['date_time'] = date_time

        out_dict['TAI93'] = []

        time93 = datetime.datetime.strptime('1993-01-01', '%Y-%m-%d')

        for i in range(len(date_time)):
            diff = date_time[i] - time93
            day_hours = 24.0
            hour_seconds = 3600.0
            out_dict['TAI93'].append(diff.days * day_hours * hour_seconds \
                    + diff.seconds)
        
        out_dict['TAI93'] = np.array(out_dict['TAI93'])

    # merge variables
    for varname in varname_list:

        # concatenate arrays
        out_dict[varname] = np.concatenate(out_dict[varname], axis=0)

        if (out_dict[varname].ndim == 4):
            # move axis
            # move lev axis to the last.
            out_dict[varname] = np.moveaxis(out_dict[varname], 1, 3)

    return out_dict 
