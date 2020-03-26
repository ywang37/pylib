"""
Created on March 24, 2020

@author: Yi Wang
"""

import datetime
import numpy as np

from mylib.io import read_nc

#
#------------------------------------------------------------------------------
#
def get_filename(root_dir, date, filename_sample):
    """
    (ywang, 03/24/20)
    Parameters
    ----------
    date : str
       'YYYY-MM-DD'
    filename_sample : str
       'GEOSFP.YYYYMMDD.A1.2x25.nc' for example

    """

    if root_dir[-1] != '/':
        root_dir = root_dir + '/'

    year = date[0:4]
    month = date[5:7]
    day = date[8:10]
    YYYYMMDD = year + month + day

    filename = root_dir + year + '/' + month + '/' + \
            filename_sample.replace('YYYYMMDD', YYYYMMDD)

    return filename
#
#------------------------------------------------------------------------------
#
def get_geosfp_hourly_A1(root_dir, date, varnames,
        filename='GEOSFP.YYYYMMDD.A1.2x25.nc', verbose=True):
    """ GEOSFP.YYYYMMDD.A1.*.nc times are 00:30,
    01:30, ......, 23:30. The function read data
    and interpolate them to 01:00, 02:00, ......,
    24:00.
    (ywang, 03/24/20)

    Parameters
    ----------
    root_dir : str
        Data root directory
    date : str
        YYYYMMDD
    varnames : list
        Variables to be processed
    fillename : str
        For example, 'GEOSFP.20180630.A1.2x25.nc'
    verbose : bool
        Output more information.

    Returns
    -------
    out_dict : dict

    """

    # 24 hourl a day
    N = 24

    if verbose:
        print(' - get_geosfp_hourly_A1: process ' + date)

    # Date
    currDate_D = datetime.datetime.strptime(date, '%Y%m%d')
    nextDate_D = currDate_D + datetime.timedelta(days=1)

    # directory
    if root_dir[-1] != '/':
        root_dir = root_dir + '/'

    # files
    curr_file = get_filename(root_dir, str(currDate_D), filename)
    next_file = get_filename(root_dir, str(nextDate_D), filename)

    # read and process data
    curr_dict = read_nc(curr_file, varnames, verbose=True)
    next_dict = read_nc(next_file, varnames, verbose=True)
    out_dict = {}
    for varn in varnames:
        tmp = np.vstack((curr_dict[varn], next_dict[varn]))
        out_dict[varn] = (tmp[0:N] + tmp[1:N+1]) * 0.5

    return out_dict
#
#------------------------------------------------------------------------------
#
def get_geosfp_hourly_A1_3days(root_dir, date, varnames,
        filename='GEOSFP.YYYYMMDD.A1.2x25.nc', 
        no_pre=False, no_next=False,
        verbose=True):
    """ Call get_geosfp_hourly_A1 to get 3-day data.
    (ywang, 03/25/20)

    Parameters
    ----------
    root_dir : str 
        Data root directory
    date : str 
        YYYYMMDD
    varnames : list
        Variables to be processed
    fillename : str
        For example, 'GEOSFP.20180630.A1.2x25.nc'
    no_pre : bool
        If true, data of the previous day is not included.
    no_next : bool
        If true, data of the next day is not included.
    verbose : bool
        Output more information.

    Returns
    -------
        out_dict : dict
    
    """

    # Date
    currDate_D = datetime.datetime.strptime(date, '%Y%m%d')
    preDate_D  = currDate_D + datetime.timedelta(days=-1)
    nextDate_D = currDate_D + datetime.timedelta(days=1)
    preDate  = str(preDate_D)
    nextDate = str(nextDate_D)

    p_date = preDate[0:4]  + preDate[5:7]  + preDate[8:10]
    c_date = date
    n_date = nextDate[0:4] + nextDate[5:7] + nextDate[8:10]

    # all dates in a list
    date_list = [p_date, c_date, n_date]
    if no_pre:
        date_list.pop(0)
    if no_next:
        date_list.pop(-1)


    # get all data
    out_dict = {}

    for varn in varnames:
        out_dict[varn] = []

    for i in range(len(date_list)):

        tmp_dict = get_geosfp_hourly_A1(root_dir, date_list[i], varnames,
                filename=filename, verbose=verbose)

        for varn in varnames:
            out_dict[varn].append(tmp_dict[varn])

    for varn in varnames:
        out_dict[varn] = np.vstack(out_dict[varn])

    return out_dict














