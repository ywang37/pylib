"""
Created on January 13, 2020

@author: Yi Wang
"""

import datetime
import numpy as np
import os
import PseudoNetCDF as pnc

def read_nd49_resample(file_list, varname_list):
    """ read data from GEOS-Chem ND49 file from resampling 
    according to satellite overpass time.

    Paratemters
    -----------
    file_list: list
        A list of ND49 filenames. It usually includes files
        of previous one day, current day, and next one day. 
    varname_list : list
        A list of variable names.

    Returns
    -------

    """

    out_dict = {}

    # 4-D variables
    for varname in varname_list:
        out_dict[varname] = []

    # other variables
    out_dict['time'] = []

    # read data
    latlon_flag = True
    for i in range(len(file_list)):

        filename = file_list[i]

        # Determine if file exists
        if not os.path.exists(filename):
            print(' - read_nd49_resample: WARNING! ' + filename + \
                    ' does not exist.')
            continue

        # read file
        print('reading ' + filename)
        infile = pnc.pncopen(filename, format='bpch')

        # get variables
        for varname in varname_list:
            out_dict[varname].append(np.array(infile.variables[varname]))

        # other variables
        out_dict['time'].append(np.array(infile.variables['time']))

        if latlon_flag:
            out_dict['latitude']  = np.array(infile.variables['latitude'])
            out_dict['longitude'] = np.array(infile.variables['longitude'])
            lat_e = np.array(infile.variables['latitude_bounds'])
            lon_e = np.array(infile.variables['longitude_bounds'])
            out_dict['latitude_e']  = np.append(lat_e[:,0], lat_e[-1,1])
            out_dict['longitude_e'] = np.append(lon_e[:,0], lon_e[-1,1])

    # merge variables
    for varname in varname_list:

        # concatenate arrays
        out_dict[varname] = np.concatenate(out_dict[varname], axis=0)

        # move axis
        out_dict[varname] = np.moveaxis(out_dict[varname], 1, 3)

    # tau
    # 'time': hours since 1985-01-01T00:00:00
    # 'TAI93': seconds since 1993-01-01T00:00:00
    # 'Time' : datetime instance
    day_hours = 24.0
    hour_seconds = 3600.0
    day_seconds = day_hours * hour_seconds
    out_dict['time'] = np.hstack(out_dict['time'])
    diff_93_85 = datetime.datetime.strptime('1993-01-01', '%Y-%m-%d') \
            - datetime.datetime.strptime('1985-01-01', '%Y-%m-%d')
    diff_93_85_hours = diff_93_85.days * day_hours
    TAI93 = (out_dict['time'] - diff_93_85_hours) * hour_seconds
    out_dict['TAI93'] = TAI93
    out_dict['Time'] = []
    for i in range(len(out_dict['TAI93'])):
        out_dict['Time'].append( \
                datetime.datetime.strptime('1993-01-01', '%Y-%m-%d') \
                + datetime.timedelta(seconds=out_dict['TAI93'][i]) )

    return out_dict

def read_bpch(filename, varname_list, latlon_flag=True,
        verbose=True, squeeze=False):
    """ read data from GEOS-Chem bpch file

    Paratemters
    -----------
    filename : str
        bpch file
    varname_list : list
        A list of variable names.
    latlon_flag : bool
        Whether latitude and longitude is returned
    squeeze : logical
        Remove single-dimensional entries from the shape
        of an array.

    Returns
    -------

    """

    out_dict = {}



    # read file
    if verbose:
        print('reading ' + filename)
    infile = pnc.pncopen(filename, format='bpch')

    # get variables
    for varname in varname_list:
        out_dict[varname] = np.array(infile.variables[varname])

    if latlon_flag:
        out_dict['latitude']  = np.array(infile.variables['latitude'])
        out_dict['longitude'] = np.array(infile.variables['longitude'])
        lat_e = np.array(infile.variables['latitude_bounds'])
        lon_e = np.array(infile.variables['longitude_bounds'])
        out_dict['latitude_e']  = np.append(lat_e[:,0], lat_e[-1,1])
        out_dict['longitude_e'] = np.append(lon_e[:,0], lon_e[-1,1])

    # Remove single-dimensional entries from the shape
    # of an array
    if squeeze:
        for varn in out_dict:
            out_dict[varn] = np.squeeze(out_dict[varn])

    return out_dict



