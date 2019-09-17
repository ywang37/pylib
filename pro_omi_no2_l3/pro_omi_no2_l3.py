"""
Created on September 15, 2019

@author: Yi Wang
"""

import glob
import numpy as np

from mylib.read_omi import read_OMI_NO2_L3

from mylib.pro_omi_no2_l3.io_omi_no2_l3 \
        import output_multi_year_month_OMI_NO2_L3
from mylib.pro_omi_no2_l3.io_omi_no2_l3 import read_month_OMI_NO2_L3



def cal_month_OMI_NO2_L3(yyyymm, in_dir, verbose=False):
    """ Calculate monthly average of OMI L3 NO2.

    Parameters:
    yyyymm : str
        Year and month. For example, 201908
    in_dir: str
        Directory for OMI L3 NO2 files
    verbose : logical
        Output more information

    Returns:
    None or ave_NO2 : None or 2D numpy array
        ave_NO2 is monthly mean NO2

    """
    min_valid = -1e30

    yyyy = yyyymm[0:4]
    mm   = yyyymm[4:6]

    # get all filenames of this month
    wildcard = in_dir + 'OMI-Aura_L3-OMNO2d_' + yyyy + 'm' + mm \
            + '*.he5'
    all_files = glob.glob(wildcard)
    all_files.sort()

    if len(all_files) == 0:
        print(' - cal_month_omi_no2_l3: WARNING, there is no file in ' 
                + yyyymm)
        return None

    # read all data
    all_NO2 = []
    for filename in all_files:

        data = read_OMI_NO2_L3(filename, verbose=verbose)
        NO2 = data['ColumnAmountNO2TropCloudScreened']
        all_NO2.append(NO2)

    all_NO2 = np.array(all_NO2)
    all_NO2[all_NO2<min_valid] = np.nan

    # monthly mean NO2
    ave_NO2 = np.nanmean(all_NO2, axis=0)

    return ave_NO2

def cal_multi_year_month_OMI_NO2_L3(month, 
        start_year, end_year, in_dir, out_dir=None, verbose=False):
    """ Calculate multi-year mean of monthly  OMI L3 NO2.

    Parameters:
    month : str
        For example, '01', '02', ...
    start_year : int
    end_year : int
    in_dir : str
    out_dir : str or None
        If it is str, output data in out_dir.
        If it is None, don't output data.
    verbose : logical
        Output more information
         
    Returns:
    ave_NO2 : 2D numpy array
        multi-year mean of monthly  OMI L3 NO2

    """

    all_NO2 = []
    for year in range(start_year, end_year+1):

        year_c = str(year)
        in_filename = in_dir + 'OMI_NO2_' + year_c + '-' +  month \
                + '_monthly.nc' 

        month_NO2 = read_month_OMI_NO2_L3(in_filename, verbose=verbose)
        all_NO2.append(month_NO2)

    all_NO2 = np.array(all_NO2)

    # mean
    ave_NO2 = np.nanmean(all_NO2, axis=0)

    # output data
    if out_dir is not None:
        output_multi_year_month_OMI_NO2_L3(out_dir, 
                str(start_year), str(end_year), month, 
                ave_NO2, verbose=verbose)

    return ave_NO2 
