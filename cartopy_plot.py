"""
Created on August 29, 2019

@author: Yi Wang
"""

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import copy
import matplotlib.pyplot as plt
import numpy as np

from mylib.pro_satellite import scale_image_multi_c

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
        cmap=plt.get_cmap('rainbow'), bad_c='grey', bad_a=1.0, 
        title=None, enhance=True, **kwargs):
    """ (1) Transfer some default parameters to pcolormesh.
        (2) Let pcolormesh can plot RGB image when *C* is a
            3D array.

    Parameters
    ----------
    ax : GeoAxes
    X : 
        Longitude
    Y : 
        Latitude
    C : 2D or 3D numpy array
        2D: (lat_dim, lon_dim)
        3D: (lat_dim, lon_dim, 3), the third dim is
            for R, G, and B channels.

    valid_min :
        Values smaller than valid_min is masked.
    valid_max :
        Values larger than valid_max is masked.
    cmap :
    bad_c : str
        The color for masked pixels
    bad_a : float
        Transparency of masked pixels
    title :
    enhance : logical
        Image enhancement

    Returns
    -------
    mesh :

    """

    # pseudocolor plot
    if (C.ndim == 2):

        # Mask array
        C_ma = copy.deepcopy(C)
        if valid_min is not None:
            C_ma = np.ma.masked_array(C_ma, C_ma<valid_min)
        if valid_max is not None:
            C_ma = np.ma.masked_array(C_ma, C_ma>valid_max)

        # The color that represents masked values
        cmap.set_bad(bad_c, alpha=bad_a)
   
        mesh = ax.pcolormesh(X, Y, C_ma, cmap=cmap, transform=ax.projection,
                **kwargs)

    # true color image
    else:

        mesh_rgb = C[:, :-1, :]

        if enhance:
            mesh_rgb = scale_image_multi_c(mesh_rgb)

        colorTuple = \
                mesh_rgb.reshape((mesh_rgb.shape[0] * mesh_rgb.shape[1]), 3)
        colorTuple = np.insert(colorTuple,3,1.0,axis=1)
        mesh = ax.pcolormesh(X, Y, C[:,:,0], color=colorTuple, 
                transform=ax.projection)

    # title
    if title is not None:
        ax.set_title(title)

    return mesh

def scatter(ax, X, Y, **kwargs):
    """ Transfer some default parameters to scatter.

    Parameters
    ----------
    ax : GeoAxes
    X :
        Longitude
    Y :
        Latitude

    """

    ax.scatter(X, Y, **kwargs)

def plot_polygon(ax, corner_lat, corner_lon, title=None, **kwargs):
    """ Plot satellite granule contour

    Parameters
    ----------
    ax : GeoAxes
    corner_lat : array-like, shape (n_points,)
        Corner latitudes
    corner_lon : array-like, shape (n_points,)
        Corner longitudes
    title : str
        Title label

    """

    n_points = corner_lat.shape[0]

    corner_lat_c = np.zeros((n_points+1,))
    corner_lat_c[0:n_points] = corner_lat
    corner_lat_c[n_points] = corner_lat[0]

    corner_lon_c = np.zeros((n_points+1,))
    corner_lon_c[0:n_points] = corner_lon
    corner_lon_c[n_points] = corner_lon[0]

    ax.plot(corner_lon_c, corner_lat_c, transform=ccrs.Geodetic(), **kwargs)
