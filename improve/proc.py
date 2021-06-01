"""
Created on May 17, 2021

@author: Yi Wang
"""

import numpy as np
import pandas as pd

#
#------------------------------------------------------------------------------
#
def split_improve(df):
    """ Split original IMPROVE data by sites.
    (ywang, 05/17/2021)

    Parameters
    ----------
    df : DataFrame
        Original IMPROVE data

    Returns
    -------
    out_dict : dict
        Keys are SiteCode
        Values are dataFrame

    """

    # Extract site information
    SiteCode = df['SiteCode']
    SiteCode = set(SiteCode)
    SiteCode = list(SiteCode)
    SiteCode.sort()

    # Split data
    out_dict = {}
    for one_site in SiteCode:
        one_df = df.loc[df['SiteCode'] == one_site]
        one_df = one_df.reset_index(drop=True)
        out_dict[one_site] = one_df

    return out_dict
#
#------------------------------------------------------------------------------
#
def extct_improve(df, cols_in=[]):
    """ Extract columns.
    (ywang, 05/17/2021)

    Parameters
    ----------
    df : DataFrame
        IMPROVE data.
    cols_in : list
        columns to be extracted.

    Returns
    -------
    out_df : DataFrame
        Extracted data.

    """

    cols = ['SiteCode', 'Date', 'Latitude', 'Longitude', 'Elevation']

    cols.extend(cols_in)

    out_df = df[cols]

    return out_df
#
#------------------------------------------------------------------------------
#
def calc_monthly_ave(date_arr, val_arr, 
        start_year=1988, end_year=2020, num_thre=3,
        undef=-999.0):
    """ Calculate monthly average.
    (ywang, 05/18/2021)

    Parameters
    ----------
    date_arr : 1-D numpy int array
        Elements are yyyymmdd
    val_arr : 1-D numpy float array
        Species concentrations
    start_year : int
        Start year
    end_year : int
        End year
    num_thre : int
        Number of days threshold for a month
    undef : float
        Used if monthly average is not available

    Returns
    -------
    out_dict : dict
        keys are 'month' and 'ave'

    """

    st_mon  = 1
    end_mon = 12

    year_arr = []
    mon_arr  = []
    day_arr  = []
    for i in range(len(date_arr)):

        date_int = date_arr[i]

        date_c = str(date_int)
        year_arr.append(date_c[0:4])
        mon_arr.append(date_c[4:6])
        day_arr.append(date_c[6:8])

    year_arr = np.array(year_arr)
    mon_arr  = np.array(mon_arr)
    day_arr  = np.array(day_arr)

    out_month = []
    out_ave   = []

    for iyr in range(start_year, end_year+1):

        curr_year = str(iyr)

        for imo in range(st_mon, end_mon+1):

            curr_mon = str(imo).zfill(2)

            out_month.append(curr_year + curr_mon)

            # get data in a month
            flag = np.logical_and(year_arr == curr_year,
                    mon_arr == curr_mon)
            flag = np.logical_and(flag,
                    val_arr > 0.0)

            # monthly ave
            num = np.sum(flag)
            if (num >= num_thre):
                tmp_data = val_arr[flag]
                ave = np.mean(tmp_data)
            else:
                ave = undef
            out_ave.append(ave)

    out_month = np.array(out_month)
    out_ave   = np.array(out_ave)

    out_dict = {}
    out_dict['month'] = out_month
    out_dict['ave']   = out_ave

    return out_dict
#
#------------------------------------------------------------------------------
#
