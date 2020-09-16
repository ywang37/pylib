"""
Created on August 29, 2019

@author: Yi Wang
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import copy
import matplotlib.pyplot as plt
import numpy as np

from pro_satellite import calculate_pixel_edge2
from pro_satellite import scale_image_multi_c

#
#------------------------------------------------------------------------------
#

def add_geoaxes(fig, *args,
        xtick=None, ytick=None, zero_direction_label=False,
        dateline_direction_label=False, number_format='g',
        degree_symbol=u'\u00B0', 
        cl_res='110m',
        cl_color=None,
        lw=None,
        title=None,
        countries=False,
        states=False,
        **kwargs):
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
    cl_res : str
        Coastline resolution. 
        Currently can be one of “110m”, “50m”, and “10m”
    title : str or None(default)
        title
    countries : bool
        Plot countries
    states : bool
        Plot states and provinces

    Returns
    -------
    ax : GeoAxes

    """

    # Set some default variables
    if xtick is None:
        xtick = np.arange(-180, 180.1, 60)
    if ytick is None:
        ytick = np.arange(-90, 90.1, 30)

    # Default projection
    kwargs['projection'] = kwargs.get('projection', ccrs.PlateCarree())
    crs = kwargs.get('projection')

    ax = fig.add_subplot(*args, **kwargs)

    if cl_color is None:
        cl_color = 'black'
    ax.coastlines(resolution=cl_res, color=cl_color, lw=lw)

    if countries:
        ax.add_feature(cfeature.BORDERS)

    if states:
        states_provinces = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='50m',
                facecolor='none')
        ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)

    # Tick labels
    tick_proj = ['PlateCarree', 'Mercator']
    proj = str(type(ax.projection)).split('.')[-1][0:-2]
    if proj in tick_proj:
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

    # title
    if title is not None:
        ax.set_title(title)

    return ax

#
#------------------------------------------------------------------------------
#

def pcolormesh(ax, X, Y, C, valid_min=None, valid_max=None, 
        cmap=None, bad_c='darkgrey', bad_a=None, 
        enhance=True, cbar=False, cbar_prop=dict(),
        **kwargs):
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
    enhance : logical
        Image enhancement
    cbar : logical
        Plot colorbar
    cbar_prop : dict
        Proterty for plt.colorbar

    Returns
    -------
    out_dict : dict

    """

    out_dict = dict()

    # pseudocolor plot
    if (C.ndim == 2):

        # Mask array
        C_ma = copy.deepcopy(C)
        if valid_min is not None:
            C_ma = np.ma.masked_array(C_ma, C_ma<valid_min)
        if valid_max is not None:
            C_ma = np.ma.masked_array(C_ma, C_ma>valid_max)

        if cmap is None:
            cmap = plt.get_cmap('rainbow')
        # The color that represents masked values
        if bad_c is not None:
            cmap.set_bad(color=bad_c, alpha=bad_a)
   
        mesh = ax.pcolormesh(X, Y, C_ma, cmap=copy.deepcopy(cmap), 
                transform=ccrs.PlateCarree(),
                **kwargs)
        out_dict['mesh'] = mesh

        # colorbar
        ax.reset_position()
        if cbar:
            cb = plt.colorbar(mesh, ax=ax, **cbar_prop)
            out_dict['cb'] = cb

    # true color image
    else:

        mesh_rgb = C[:, :-1, :]

        if enhance:
            mesh_rgb = scale_image_multi_c(mesh_rgb)

        colorTuple = \
                mesh_rgb.reshape((mesh_rgb.shape[0] * mesh_rgb.shape[1]), 3)
        colorTuple = np.insert(colorTuple,3,1.0,axis=1)
        mesh = ax.pcolormesh(X, Y, C[:,:,0], color=colorTuple, 
                transform=ccrs.PlateCarree())
        out_dict['mesh'] = mesh


    return out_dict

#
#------------------------------------------------------------------------------
#
