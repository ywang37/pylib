"""
Created on March 19, 2020

@author: Yi Wang
"""

from copy import deepcopy
from netCDF4 import Dataset
import numpy as np
import os

from mylib.io import read_nc

#
#------------------------------------------------------------------------------
#
def download_met(res, model, year, month, root_dir=None):
    """ Download meteorological data for GEOS-Chem
    (ywang, 03/19/2020)

    Parameters
    ----------
    res : str
        Resolution. 'GEOS_2x2.5', 'GEOS_0.25x0.3125_CH', for example.
    model : str
        'GEOS_FP', 'MERRA2'
    year : str
        '2018' for example.
    month : str
        '05' for example
    root_dir : str or None
        Parent directory to save data.
        If None, root_dir is '/Dedicated/jwang-data/GCDATA/'

    """

    # root URL
    url_root = 'http://geoschemdata.computecanada.ca/ExtData/'

    # URL
    url = url_root + res + '/' + model + '/' + year + '/' + \
            month + '/'
    #print(url)

    # root directory
    if root_dir is None:
        root_dir = '/Dedicated/jwang-data/GCDATA/'
    if root_dir[-1] != '/':
        root_dir = root_dir + '/'

    # data directory
    data_year_dir = root_dir + res + '/' + model + '/' + year + '/' 
    data_month_dir = data_year_dir + month + '/'

    # check year directory
    if not os.path.isdir(data_year_dir):
        os.system('mkdir -p ' + data_year_dir)
        os.system('chmod g+w ' + data_year_dir)

    # check month directory
    if not os.path.isdir(data_month_dir):
        os.system('mkdir -p ' + data_month_dir)
        os.system('chmod g+w ' + data_month_dir)
    else:
        print(data_month_dir + ' already exists, data is not downloaded.')
        return

    # download
    cmd = 'wget -r -np -nH -N -R "*.html*" --cut-dirs=5 ' + \
            '"' + url + '" -P "' + data_month_dir + '"'
    os.system(cmd)

    # access right
    os.system('chmod g+w ' + data_month_dir + '*')

#
#------------------------------------------------------------------------------
#
def link_for_soil_temp(res, model, year, month, root_dir=None):
    """ All meteorological data are lined except 'A1' data.
    New 'A1' data contain soil temperatuer.
    (ywang, 03/19/2020)

    Parameters
    ----------
    res : str
        Resolution. 'GEOS_2x2.5', 'GEOS_0.25x0.3125_CH', for example.
    model : str
        'GEOS_FP', 'MERRA2'
    year : str
        '2018' for example.
    month : str
        '05' for example
    root_dir : str or None
        Parent directory for meterorological data.
        If None, root_dir is '/Dedicated/jwang-data/GCDATA/'

    """

    if model == 'GEOS_FP':
        model_file = 'GEOSFP'
        nc_suffix = 'nc'
    elif model == 'MERRA2':
        model_file = 'MERRA2'
        nc_suffix = 'nc4'

    # root directory
    if root_dir is None:
        root_dir = '/Dedicated/jwang-data/GCDATA/'
    if root_dir[-1] != '/':
        root_dir = root_dir + '/'

    # Original data directory
    ori_data_dir = root_dir + res + '/' + model + '/' + year + '/' \
            + month + '/'
    print('Original directory is ' + ori_data_dir)

    # New data directory
    new_data_dir = root_dir + res + '/' + model + '_soil_T/' + year + '/' \
            + month + '/'
    print('New directory is ' + new_data_dir)

    # check directory
    if not (os.path.isdir(ori_data_dir)  and os.path.isdir(new_data_dir)):
        print('At least one directory does not exist.')
        print('Link is not done.')

    # link
    os.system('ln -sf ' + ori_data_dir + model_file + '.' + year + month +
            '*.A3*.*.' + nc_suffix  + ' ' + new_data_dir)
    os.system('ln -sf ' + ori_data_dir + model_file + '.' + year + month +
            '*.I3.*.' + nc_suffix + ' ' + new_data_dir)
#
#------------------------------------------------------------------------------
#
def correct_BC(curr_file, pre_file, old_times=7):
    """ Boundary conditions at the fisrt time is missing, we need
    to fill the gap.
    (ywang, 08/26/2020)

    Parameters
    ----------
    curr_file : str
        The file to be updated
    pre_file : str or None
        The file for boundary conditions at the fisrt time
    old_times : int
        The code execcute if the times in curr_file equals old_times

    """

    print(' - correct_BC: curr_file is: ' + curr_file)
    print(' - correct_BC: pre_file is: ', pre_file)

    # check time
    curr_tmp = read_nc(curr_file, ['time'])
    curr_times = curr_tmp['time'].shape[0]
    if curr_times != old_times:
        print(' - correct_BC: curr_times={}'.format(curr_times))
        print(' - correct_BC: old_times={}'.format(old_times))
        print(' - correct_BC: just skip.')
        return

    # copy file
    curr_file_old = curr_file.split('.')
    curr_file_old[-2] = curr_file_old[-2] + '_old'
    curr_file_old = '.'.join(curr_file_old)
    if not os.path.exists(curr_file_old):
        os.system('cp ' + curr_file + ' ' + curr_file_old)

    # open curr_file
    f_curr = Dataset(curr_file, 'r+')

    # open pre_file
    if pre_file is not None:
        f_pre = Dataset(pre_file, 'r')
        # check data
        curr_time_units = f_curr['time'].units
        pre_time_units  = f_pre['time'].units
        if (curr_time_units != curr_time_units):
            print('curr_time_units is: ' + curr_time_units)
            print('pre_time_units is: ' + pre_time_units)
            exit()

    # fill time
    print(' - correct_BC: process time')
    for i in np.array(range(curr_times))[::-1]:
        f_curr['time'][i+1] = f_curr['time'][i]
    if pre_file is not None:
        f_curr['time'][0] =  f_pre['time'][0]
    else:
        f_curr['time'][0] =  0

    # fill data
    all_varns = list(f_curr.variables.keys())
    for varn in all_varns:
        if 'SpeciesBC_' in varn:
            print(' - correct_BC: process ' + varn)
            for i in np.array(range(curr_times))[::-1]:
                f_curr[varn][i+1,:,:,:] = f_curr[varn][i,:,:,:]
            if pre_file is not None:
                f_curr[varn][0,:,:,:] =  f_pre[varn][0,:,:,:]

    # close files
    f_curr.close()
    if pre_file is not None:
        f_pre.close()
#
#------------------------------------------------------------------------------
#
