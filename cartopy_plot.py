"""
Created on August 29, 2019

@author: Yi Wang
"""

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt
import numpy as np

def add_geoaxes(fig, *args, xtick=np.arange(-180, 180.1, 60), 
        ytick=np.arange(-90, 90.1, 30), zero_direction_label=False,
        dateline_direction_label=False, number_format='g',
        degree_symbol=u'\u00B0', **kwargs):
    """ Add a GeoAxes instance to Figure (fig) instance.

    Parameters
    ----------
    fig : Figure
    *args
        It is idential to *args in fig.add_subplot
    **kwargs
        It is idential to **kwargs in fig.add_subplot
    xtick :
        Longitude
    ytick :
        Latitude
    zero_direction_label :
        Direction label at 0 degree longitude
    dateline_direction_label :
        Direction label at 180 degree longitude
    number_format :
    degree_symbol : 

    Returns
    -------
    ax : GeoAxes

    """

    # Default projection
    kwargs['projection'] = kwargs.get('projection', ccrs.PlateCarree())
    crs = kwargs.get('projection')

    ax = fig.add_subplot(*args, **kwargs)

    ax.coastlines()

    # Tick labels
    ax.set_xticks(xtick, crs=ccrs.PlateCarree())
    ax.set_yticks(ytick, crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(
            zero_direction_label=zero_direction_label,
            dateline_direction_label=dateline_direction_label,
            number_format=number_format, degree_symbol=degree_symbol)
    lat_formatter = LatitudeFormatter(number_format=number_format,
            degree_symbol=degree_symbol)
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    return ax

def pcolormesh(ax, X, Y, C, valid_min=None, valid_max=None, 
        cmap=plt.get_cmap('rainbow'), title=None, **kwargs):
    """ Transfer some default parameters to pcolormesh.

    Parameters
    ----------
    ax : GeoAxes
    X : 
        Longitude
    Y : 
        Latitude
    C : 
        The values will be colored
    valid_min :
        Values smaller than valid_min is masked.
    valid_max :
        Values larger than valid_max is masked.
    title :

    Returns
    -------
    mesh :

    """

    # Mask array
    C_ma = np.array(C)
    if valid_min is not None:
        C_ma = np.ma.masked_array(C_ma, C_ma<valid_min)
    if valid_max is not None:
        C_ma = np.ma.masked_array(C_ma, C_ma>valid_max)

    # The color that represents masked values
    cmap.set_bad('grey')

    # title
    if title is not None:
        ax.set_title(title)

    mesh = ax.pcolormesh(X, Y, C_ma, cmap=cmap, **kwargs)

    return mesh












