"""
Created on November 6, 2019

@author: Yi Wang
"""

import numpy as np

from mylib.cartopy_plot import cartopy_plot
from mylib.pro_tropomi.read_tropomi import read_tropomi_no2

#
#------------------------------------------------------------------------------
#

def plot_granule_tropomi_no2(*args, ax=None, fig=None,
        xtick=None, ytick=None,
        region_limit=None,
        valid_min=None, valid_max=None,
        cbar=True, cbar_prop={},
        title='',
        **kwargs,
        ):
    """ Plot a granule of TROPOMI
    nitrogendioxide_tropospheric_column.

    Parameters
    ----------
    Most parameters are transfered to
    cartopy_plot in cartopy_plot.py


    """

    # get *agrs
    if (len(args) == 4):
        filename = args[3]
        var_name = args[0:3]
        in_data = read_tropomi_no2(filename, varnames=var_name, 
                replace=True, verbose=True)
        lat = in_data[args[0]]
        lon = in_data[args[1]]
        var = in_data[args[2]]
    else:
        lat = copy.deepcopy(args[0])
        lon = copy.deepcopy(args[1])
        var = copy.deepcopy(args[2])

    out_dict = cartopy_plot(lat, lon, var, ax=None, fig=None,
            xtick=xtick, ytick=ytick,
            region_limit=region_limit,
            valid_min=valid_min, valid_max=valid_max,
            cbar=cbar, cbar_prop=cbar_prop,
            title=title,
            **kwargs)

    return out_dict

#
#------------------------------------------------------------------------------
#
