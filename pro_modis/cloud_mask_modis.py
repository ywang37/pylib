"""
Created on December 05, 2019

@author: Yi Wang
"""

import numpy as np

#
#------------------------------------------------------------------------------
#
def stddev_3x3(RefSb, flag_cal=None):
    """ Calculate 3x3 box standard deviation.
    (Yi Wang, 12/05/2019)

    Parameters
    ----------
    RefSb : 2-D numpy array
        TOA reflectance

    flag_cal : 2-D numpy array
        Only process pixels with flag_cal is True.
        If flag_cal is None, all elements are set as True

    """

    undef = -999.0

    # 3x3 box
    step = 3

    # standard deviation array
    stddev_a = np.full_like(RefSb, undef)

    # average array
    ave_a = np.full_like(RefSb, undef)

    # False means standard deviation is not calculated.
    flag = np.full_like(RefSb, False)

#    # flag_cal
#    if flag_cal is None:


    #---------------------------------
    # 3x3 box standard deviation
    #---------------------------------
    dim = RefSb.shape
    nline  = dim[0]
    npixel = dim[1]
    for i in range(0, nline, step):
        for j in range(0, npixel, step):

            # index for a 3x3 box
            in_i1 = i
            in_i2 = in_i1 + step - 1
            in_j1 = j
            in_j2 = in_j1 + step - 1

            # standard deviation
            stddev = np.std(RefSb[in_i1:in_i2+1, in_j1:in_j2+1])
            stddev_a[in_i1:in_i2+1, in_j1:in_j2+1] = stddev

            # standard deviation
            ave = np.mean(RefSb[in_i1:in_i2+1, in_j1:in_j2+1])
            ave_a[in_i1:in_i2+1, in_j1:in_j2+1] = ave


            # flag
            if not np.isnan(stddev):
                flag[in_i1:in_i2+1, in_j1:in_j2+1] = True
            else:
                print(stddev)

    out_dict = {}
    out_dict['stddev'] = stddev_a
    out_dict['ave']    = ave_a
    out_dict['flag']   = flag

    return out_dict
#
#------------------------------------------------------------------------------
#
def cloud_ocean_stddev(RefSb_466nm, RefSb_554nm, RefSb_645nm,
        thre_stddev=0.0025, thre_dust=0.75):
    """ 3x3 box standard deviation test over ocean.
    (Yi Wang, 12/05/2019)

    """

    # 1 represents cloud, and 0 represents clear
    cloud_mask = np.full_like(RefSb_554nm, 1, np.int)

    # RefSb_554nm standard deviation test
    stddev_dict = stddev_3x3(RefSb_554nm)
    stddev = stddev_dict['stddev']
    stddev_flag = stddev_dict['flag']
    # 0 represents clear
    cloud_mask[np.logical_and(stddev_flag, stddev<thre_stddev)] = 0

    # dust call back
    ratio = RefSb_466nm / RefSb_645nm
    print('cloud_ocean_stddev 01')
    cloud_mask[ratio<thre_dust] = 0
    print('cloud_ocean_stddev 02')

    # output
    out_dict = {}
    out_dict['cloud_mask']  = cloud_mask
    out_dict['ratio']       = ratio
    out_dict['stddev']      = stddev
    out_dict['stddev_flag'] = stddev_flag
    out_dict['ave']         = stddev_dict['ave']

    
    return out_dict
#
#------------------------------------------------------------------------------
#
