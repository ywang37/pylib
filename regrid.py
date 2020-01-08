"""
Created on November 6, 2019

@author: Yi Wang
"""

from netCDF4 import Dataset
import numpy as np
from tqdm import tqdm

from mylib.grid_utility import get_center_index_2

def drop_in_the_box(data_list, lat_list, lon_list, 
        lat_start, lat_end, lat_step,
        lon_start, lon_end, lon_step,
        weight_list=None
        ):
    """ Regrid data through the drop in the box approach.

    Parameters
    ----------
    data_list : list
        Elements are numpy array
    lat_list : list
        Elements are numpy array 
    lon_list : list
        Elements are numpy array
    lat_start : float
        Starting latitude
    lat_end : float
        Ending latitude
    lat_step : float
        Latitude step
    lon_start : float
        Starting longiitude
    lon_end : float
        Ending longitude
    lon_step : float
        Longitude step
    weight_list : list or None (defult)

    Returns
    -------
    out_data : dict
        A dictionary of all variables.
    """

    # save new data
    nlat = round( (lat_end - lat_start) / lat_step )
    nlon = round( (lon_end - lon_start) / lon_step )
    new_sum = np.zeros( (nlat, nlon) )
    new_weight = np.zeros( (nlat, nlon) )

    for d_ind in range(len(data_list)):

        print(' - drop_in_the_box: processing {}'.format(d_ind))

        # data, latitue and longitude
        data = data_list[d_ind].flatten()
        lat  = lat_list[d_ind].flatten()
        lon  = lon_list[d_ind].flatten()

        # weights
        if weight_list is None:
            weight = np.ones_like(data, dtype=float)
        else:
            weight = weight_list[d_ind].flatten()

        # subset data
        flag = (lat >= lat_start) & (lat <= lat_end) \
                & (lon >= lon_start) & (lon <= lon_end)
        data   = data[flag]
        lat    = lat[flag]
        lon    = lon[flag]
        weight = weight[flag]

        for k in tqdm(range(len(data))):

            # latitude index for the new grid
            i = get_center_index_2(lat[k], lat_start, lat_end, lat_step)

            # longitude index for the new grid
            j = get_center_index_2(lon[k], lon_start, lon_end, lon_step)

            new_sum[i,j]    += (data[k] * weight[k])
            new_weight[i,j] += weight[k]

    new_sum = np.ma.masked_array(new_sum, new_weight<=0.0)
    new_ave = new_sum / new_weight
    out_data = {}
    out_data['value'] = new_ave
    out_data['count'] = new_weight

    return out_data
