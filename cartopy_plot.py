"""
Created on August 29, 2019

@author: Yi Wang
"""

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import copy
import matplotlib.pyplot as plt
import numpy as np

from mylib.io import read_nc
from mylib.pro_satellite.pro_satellite import calculate_pixel_edge2
from mylib.pro_satellite.pro_satellite import scale_image_multi_c

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

def contourf(ax, *args, valid_min=None, valid_max=None, 
        cmap=plt.get_cmap('rainbow'), bad_c='grey', bad_a=1.0, 
        cbar=False, **kwargs):
    """ Transfer some default parameters to contourf.

    Parameters
    ----------
    ax : GeoAxes
    *args : 
        Same as *args in plt.contourf
        contourf(ax, [X, Y], Z, [levels], ...... )
    valid_min :
        Values smaller than valid_min is masked.
    valid_max :
        Values larger than valid_max is masked.
    cmap :
    bad_c : str
        The color for masked pixels
    bad_a : float
        Transparency of masked pixels
    cbar : logical
        Plot colorbar

    Returns
    -------
    out_dict : dict
        qcs: return of ax.contourf

    """

    out_dict = dict()

    # get *args
    args = list(args)
    if (len(args) == 1) or (len(args) == 2):
        C_ma = copy.deepcopy(args[0])
    elif (len(args) == 3) or (len(args) == 4):
        C_ma = copy.deepcopy(args[2])
    else:
        print(' - contourf: Number of parameters is incorrect.')
        exit()

    # mask array
    if valid_min is not None:
        C_ma = np.ma.masked_array(C_ma, C_ma<valid_min)
    if valid_max is not None:
        C_ma = np.ma.masked_array(C_ma, C_ma>valid_max)

    # put back masked array
    if (len(args) == 1) or (len(args) == 2):
        args[0] = C_ma
    elif (len(args) == 3) or (len(args) == 4):
        args[2] = C_ma

    # The color that represents masked values
    cmap.set_bad(bad_c, alpha=bad_a)
   
    qcs = ax.contourf(*args, cmap=copy.deepcopy(cmap), 
            transform=ccrs.PlateCarree(),
            **kwargs)
    out_dict['qcs'] = qcs


    # colorbar
    if cbar:
        cb = plt.colorbar(qcs, ax=ax)
        out_dict['cb'] = cb

    return out_dict

#
#------------------------------------------------------------------------------
#

def scatter(ax, X, Y,
        cbar=False, cbar_prop=dict(),
        **kwargs):
    """ Transfer some default parameters to scatter.

    Parameters
    ----------
    ax : GeoAxes
    X :
        Longitude
    Y :
        Latitude
    cbar : logical
        Plot colorbar
    cbar_prop : dict
        Proterty for plt.colorbar

    Returns
    -------
    out_dict : dict
        keys : paths
        operatioanl keys: cb

    """

    out_dict = dict()

    # scatter plot
    paths = ax.scatter(X, Y, transform=ccrs.PlateCarree(), **kwargs)
    out_dict['paths'] = paths

    # colorbar
    if cbar:
        cb = plt.colorbar(paths, ax=ax, **cbar_prop)
        out_dict['cb'] = cb

    return out_dict

#
#------------------------------------------------------------------------------
#

def plot_polygon(ax, corner_lat, corner_lon, **kwargs):
    """ Plot satellite granule contour

    Parameters
    ----------
    ax : GeoAxes
    corner_lat : array-like, shape (n_points,)
        Corner latitudes
    corner_lon : array-like, shape (n_points,)
        Corner longitudes

    """

    n_points = corner_lat.shape[0]

    corner_lat_c = np.zeros((n_points+1,))
    corner_lat_c[0:n_points] = corner_lat
    corner_lat_c[n_points] = corner_lat[0]

    corner_lon_c = np.zeros((n_points+1,))
    corner_lon_c[0:n_points] = corner_lon
    corner_lon_c[n_points] = corner_lon[0]

    ax.plot(corner_lon_c, corner_lat_c, transform=ccrs.Geodetic(), **kwargs)

def plot_rectangle(ax, region_limit, **kwargs):
    """ Plot a rectangle accoring to *region_limit*

    Parameters
    ----------
    ax : GeoAxes
    region_limit : list-like
        [lat_min, lon_min, lat_max. lon_max]

    """

    corner_lat = np.array( [region_limit[0], region_limit[0], 
            region_limit[2], region_limit[2]] )
    corner_lon = np.array( [region_limit[1], region_limit[3], 
            region_limit[3], region_limit[1]] )

    plot_polygon(ax, corner_lat, corner_lon, **kwargs)

#
#------------------------------------------------------------------------------
#

def cartopy_plot(*args, ax=None, fig=None,
        indices=None,
        reader=None, reader_prop={},
        xtick=None,
        ytick=None,
        cl_res=None,
        region_limit=None,
        valid_min=None, valid_max=None,
        cbar=True, cbar_prop = {},
        title=None,
        **kwargs):
    """ Plot a 2-D variable by pcolormesh.

    Parameters
    ----------
    *agrs:
        3 elements: longitude, latitude, variable
        4 elelemts: longitude_name, latitude_name,
                    variable_name, filename
    ax : GeoAxes or None (default)
        Create a GeoAxes if ax is None. 
    fig : plt.figure() or None (default)
        Create a plt.figure() if both fig and ax are None
    indices : None
        Used to get subset if the dimension of variable is 
        larger than 2. (The function is not implemented yet.)
    reader : function to read data.
        if reader is None, read_nc is used
    reader_prop : dict
        optional variables to reader if reader is not None
    xtick : list-like
        Longitude ticks
    ytick : list-like
        Latitude ticks
    cl_res : str
        Coastline resolution. 
        Currently can be one of “110m”, “50m”, and “10m”
    region_limit : tuple-like or None
        (min_lat, min_lon, max_lat, max_lon)
    valid_min : float or None
        Values that are smaller than valid_min is masked.
    valid_max : float or None
        Values that are larger that valid_max is masked.
    cbar : logical
        Whether or not plot colorbar.
    cbar_prpo : dict
        Colorbar properties, transferred to plt.colorbar()
    title : str
        Title
    **kwargs : dict
        Keywords to pcolormesh

    Returns
    -------
    out_dict :
        keys : ax, fig, mesh
        operatioanl keys: cb

    """

    out_dict = {}

    # get *args
    if (len(args) == 4):
        filename = args[3]
        var_name = args[0:3]
        if reader is None:
            reader = read_nc
            in_data = reader(filename, varnames=var_name, verbose=True)
        else:
            read_prop['verbose'] = read_prop.get('verbose', True)
            in_data = reader(filename, varnames=var_name, **read_prop)
        lon = in_data[args[0]]
        lat = in_data[args[1]]
        var = in_data[args[2]]
        if indices is not None:
            pass
    else:
        lon = copy.deepcopy(args[0])
        lat = copy.deepcopy(args[1])
        var = copy.deepcopy(args[2])

    # coastline resolution
    if cl_res is None:
        cl_res='110m'

    # GeoAxes
    if ax is None:
        if fig is None:
            fig = plt.figure()
        ax = add_geoaxes(fig, cl_res=cl_res, xtick=xtick, ytick=ytick)

    # set region limit
    if region_limit is not None:
        ax.set_xlim(region_limit[1], region_limit[3])
        ax.set_ylim(region_limit[0], region_limit[2])

    # title
    if title is not None:
        ax.set_title(title)

    out_dict['ax'] = ax
    out_dict['fig'] = fig

    # plot true color image
    if (var.ndim == 3) and (len(args) == 3):
        pcolormesh(ax, lon, lat, var)
        return

    # calculate edge or not
    if lat.shape == var.shape:
        lat, lon = calculate_pixel_edge2(lat, lon)

    # colorbar property
    #cbar_prop['orientation'] = cbar_prop.get('orientation', 'horizontal')

    # plot
    pout = pcolormesh(ax, lon, lat, var,
            valid_min=valid_min, valid_max=valid_max,
            cbar=cbar, cbar_prop=cbar_prop,
            **kwargs)
    out_dict['mesh'] = pout['mesh']
    out_dict['cb'] = pout.get('cb', None)

    return out_dict

#
#------------------------------------------------------------------------------
#

def cartopy_plot_scatter(*args, ax=None, fig=None,
        reader=None, reader_prop={},
        xtick=None,
        ytick=None,
        cl_res='110m',
        region_limit=None,
        cbar=True, cbar_prop = {},
        title=None,
        **kwargs):
    """ Plot scatter map.
    (Yi Wang, 11/27/2019)

    Parameters
    ----------
    *agrs:
        2 elements: longitude, latitude
        3 elelemts: longitude_name, latitude_name,
                    filename
    ax : GeoAxes or None (default)
        Create a GeoAxes if ax is None. 
    fig : plt.figure() or None (default)
        Create a plt.figure() if both fig and ax are None
    reader : function to read data.
        if reader is None, read_nc is used
    reader_prop : dict
        optional variables to reader if reader is not None
    xtick : list-like
        Longitude ticks
    ytick : list-like
        Latitude ticks
    cl_res : str
        Coastline resolution. 
        Currently can be one of “110m”, “50m”, and “10m”
    region_limit : tuple-like or None
        (min_lat, min_lon, max_lat, max_lon)
    cbar : logical
        Whether or not plot colorbar.
    cbar_prpo : dict
        Colorbar properties, transferred to plt.colorbar()
    title : str
        Title
    **kwargs : dict
        Keywords to pcolormesh

    Returns
    -------
    out_dict :
        keys : ax, fig, paths
        operatioanl keys: cb

    """

    out_dict = {}

    # get *args
    if (len(args) == 3):
        filename = args[2]
        var_name = args[0:2]
        if reader is None:
            reader = read_nc
            in_data = reader(filename, varnames=var_name, verbose=True)
        else:
            read_prop['verbose'] = read_prop.get('verbose', True)
            in_data = reader(filename, varnames=var_name, **read_prop)
        lon = in_data[args[0]]
        lat = in_data[args[1]]
        if indices is not None:
            pass
    else:
        lon = copy.deepcopy(args[0])
        lat = copy.deepcopy(args[1])

    # GeoAxes
    if ax is None:
        if fig is None:
            fig = plt.figure()
        ax = add_geoaxes(fig, cl_res=cl_res, xtick=xtick, ytick=ytick)

    # set region limit
    if region_limit is not None:
        ax.set_xlim(region_limit[1], region_limit[3])
        ax.set_ylim(region_limit[0], region_limit[2])

    # title
    if title is not None:
        ax.set_title(title)

    out_dict['ax'] = ax
    out_dict['fig'] = fig

    # colorbar property
    #cbar_prop['orientation'] = cbar_prop.get('orientation', 'horizontal')

    # plot
    pout = scatter(ax, lon, lat,
            cbar=cbar, cbar_prop=cbar_prop,
            **kwargs)
    out_dict['paths'] = pout['paths']
    out_dict['cb'] = pout.get('cb', None)

    return out_dict

#
#------------------------------------------------------------------------------
#
