"""
Created on September 2, 2020

@author: Yi Wang
"""

import cartopy.feature as cfeature
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

from mylib.cartopy_plot import cartopy_plot, cartopy_plot_scatter
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.layout import multiFigure, h_1_ax, h_2_ax, panel_tick_label
from mylib.layout import right_center_label

#
#------------------------------------------------------------------------------
#
def region_limit_to_vertex(region_limit):
    """ Convert region_limit to vertex information.

    Parameters
    ----------
    region_limit : tuple-like
        (min_lat, min_lon, max_lat, max_lon)

    Returns
    -------
    lat : tuple
        Latitudes of vertexes
    lon : tuple
        Longitudes of vertexes

    """

    lat = (region_limit[0], region_limit[0], region_limit[2], region_limit[2])
    lon = (region_limit[1], region_limit[3], region_limit[3], region_limit[1])

    return lat, lon
#
#------------------------------------------------------------------------------
#
def plot_comparison(var1, var2, lat_e, lon_e, 
        figsize=(6,7), vmin=None, vmax=None,
        diff_vmin=None, diff_vmax=None,
        cb_diff_ticks=None,
        left=0.25, right=0.75,
        top=0.95, bottom=0.1,
        hspace=None, wspace=None,
        cl_color='k',
        lw=None,
        region_limit=None,
        xticks=np.arange(-180, 180.0, 60),
        yticks=np.arange(-90, 90.1, 30),
        title_dict={},
        cmap=deepcopy(WhGrYlRd_map),
        diff_cmap=deepcopy(gbcwpry_map),
        units='',
        cb_ratio=0.03,
        cb_y_off=-0.05,
        mask_ocean=False,
        ):
    """ Plot var1, var2 and var2 - var1.
    (Yi Wang, 09/02/2020)

    Parameters
    ----------
    var1 : 2-D numpy array
        The first variable.
    var2 : 2-D numpy array
        The second variable.
    lat_e : 2-D numpy array
        Latitude edges.
    lon_e : 2-D numpy array
        Longitude edge.
    figsize : 2-element tuple
        Figure size
    vmin : float or None
        Minimum to be plotted for var1 and var2.
        If vmin is None, it will be set as the valid minimum of var1 
        or var2, which is smaller.
    vmax : float or None
        Similar to vmin, but of maximum
    diff_vmin : float or None
        Minimum to be plotted for var2 - var1.
    diff_vmax : float or None
        Maximum to be plotted for var2 = var1.

    Returns
    -------

    """

    # define figure
    nrow = 3
    ncol = 1
    projPos = list(range(nrow * ncol))
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos,
            countries=True, states=True,
            coastlines=True, mask_ocean=mask_ocean)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']
 
    # get min and max for var1 and var2
    if vmin is None:
        vmin = min(np.min(var1), np.min(var2))
    if vmax is None:
        vmax = max(np.max(var1), np.max(var2))

    # plot var1 and var2
    var_list = [var1, var2]
    varname_list = ['var1', 'var2']
    for i in range(len(var_list)):

        ax = axes[i]
        var = var_list[i]
        title = title_dict.get(varname_list[i], None)
        pout = cartopy_plot(lon_e, lat_e, var, ax=ax, title=title,
                cbar=False, vmin=vmin, vmax=vmax, cmap=deepcopy(cmap))


    # plot var2 - var1
    diff = var2 - var1
    ax = axes[2]
    title = title_dict.get('diff', None)
    diff_pout = cartopy_plot(lon_e, lat_e, diff, ax=ax, title=title,
            cbar=False, vmin=diff_vmin, vmax=diff_vmax, 
            cmap=deepcopy(diff_cmap))

    # ticks
    panel_tick_label(axes, ncol, xticks=xticks, yticks=yticks)

    # set limit
    for ax in axes:
        if region_limit is not None:
            ax.set_xlim((region_limit[1],region_limit[3]))
            ax.set_ylim((region_limit[0],region_limit[2]))


    # colorbar for var1 and var2
    cax1 = h_1_ax(fig, axes[1],ratio=cb_ratio, y_off=cb_y_off)
    cb1 = plt.colorbar(pout['mesh'], cax=cax1,
            orientation='horizontal')
    right_center_label(cax1, units)

    # colorbar for diff
    cax2 = h_1_ax(fig, axes[2],ratio=cb_ratio, y_off=cb_y_off)
    cb2 = plt.colorbar(diff_pout['mesh'], cax=cax2,
            orientation='horizontal')
    if cb_diff_ticks is not None:
        cb2.set_ticks(cb_diff_ticks)
    right_center_label(cax2, units)
#
#------------------------------------------------------------------------------
#
def plot_comparison_2(var1, var2, lat_e, lon_e, 
        figsize=(8,6), vmin=None, vmax=None,
        diff_vmin=None, diff_vmax=None,
        rel_pct_vmin=None, rel_pct_vmax=None,
        cb_diff_ticks=None,
        cb_rel_pct_ticks=None,
        extend='neither',
        diff_extend='neither',
        rel_pct_extend='neither',
        left=0.08, right=0.80,
        top=0.95, bottom=0.1,
        hspace=None, wspace=None,
        cl_color='k',
        lw=None,
        region_limit=None,
        xticks=np.arange(-180, 180.0, 60),
        yticks=np.arange(-90, 90.1, 30),
        title_dict={},
        cmap=deepcopy(WhGrYlRd_map),
        diff_cmap=deepcopy(gbcwpry_map),
        rel_pct_cmap=deepcopy(gbcwpry_map),
        units='',
        diff_units='',
        cb_ratio=0.06,
        cb_y_off_t=-0.04,
        cb_y_off_b=-0.07,
        mask_ocean=False,
        ):
    """ Plot var1, var2 and var2 - var1.
    (Yi Wang, 10/08/2020)

    Parameters
    ----------
    var1 : 2-D numpy array
        The first variable.
    var2 : 2-D numpy array
        The second variable.
    lat_e : 2-D numpy array
        Latitude edges.
    lon_e : 2-D numpy array
        Longitude edge.
    figsize : 2-element tuple
        Figure size
    vmin : float or None
        Minimum to be plotted for var1 and var2.
        If vmin is None, it will be set as the valid minimum of var1 
        or var2, which is smaller.
    vmax : float or None
        Similar to vmin, but of maximum
    diff_vmin : float or None
        Minimum to be plotted for var2 - var1.
    diff_vmax : float or None
        Maximum to be plotted for var2 = var1.
    rel_pct_vmin : float or None
        Minimum to be plotted for (var2 - var1) / var1 * 100
    rel_pct_vmax : float or None
        Maximum to be plotted for (var2 - var1) / var1 * 100

    Returns
    -------

    """

    # define figure
    nrow = 2
    ncol = 2
    projPos = list(range(nrow * ncol))
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos,
            countries=True, states=True,
            coastlines=True, mask_ocean=mask_ocean)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']
 
    # get min and max for var1 and var2
    if vmin is None:
        vmin = min(np.nanmin(var1), np.nanmin(var2))
    if vmax is None:
        vmax = max(np.nanmax(var1), np.nanmax(var2))

    # plot var1 and var2
    var_list = [var1, var2]
    varname_list = ['var1', 'var2']
    for i in range(len(var_list)):

        ax = axes[i]
        var = var_list[i]
        title = title_dict.get(varname_list[i], None)
        pout = cartopy_plot(lon_e, lat_e, var, ax=ax, title=title,
                cbar=False, vmin=vmin, vmax=vmax, cmap=deepcopy(cmap))

    # plot (var2 - var1) / var1 * 100
    rel_pct = (var2 - var1) / var1 * 100
    if (rel_pct_vmin is None) and (rel_pct_vmax is None):
        rel_pct_vmax = np.nanmax(np.absolute(rel_pct))
        rel_pct_vmin = -rel_pct_vmax
    ax = axes[2]
    title = title_dict.get('relative_diff', None)
    rel_pct_pout = cartopy_plot(lon_e, lat_e, rel_pct, ax=ax, title=title,
            cbar=False, vmin=rel_pct_vmin, vmax=rel_pct_vmax,
            cmap=deepcopy(rel_pct_cmap))

    # plot var2 - var1
    diff = var2 - var1
    if (diff_vmin is None) and (diff_vmax is None):
        diff_vmax = np.nanmax(np.absolute(diff))
        diff_vmin = -diff_vmax
    ax = axes[3]
    title = title_dict.get('diff', None)
    diff_pout = cartopy_plot(lon_e, lat_e, diff, ax=ax, title=title,
            cbar=False, vmin=diff_vmin, vmax=diff_vmax, 
            cmap=deepcopy(diff_cmap))

    # ticks
    panel_tick_label(axes, ncol, xticks=xticks, yticks=yticks)

    # set limit
    for ax in axes:
        if region_limit is not None:
            ax.set_xlim((region_limit[1],region_limit[3]))
            ax.set_ylim((region_limit[0],region_limit[2]))

    # colorbar for var1 and var2
    cax1 = h_2_ax(fig, axes[0], axes[1],ratio=cb_ratio, y_off=cb_y_off_t)
    cb1 = plt.colorbar(pout['mesh'], cax=cax1,
            orientation='horizontal', extend=extend)
    right_center_label(cax1, units)

    # colorbar for (var2 - var1) / var1 * 100
    cax2 = h_1_ax(fig, axes[2],ratio=cb_ratio, y_off=cb_y_off_b)
    cb2 = plt.colorbar(rel_pct_pout['mesh'], cax=cax2,
            orientation='horizontal', extend=rel_pct_extend)
    if cb_rel_pct_ticks is not None:
        cb2.set_ticks(cb_rel_pct_ticks)
    right_center_label(cax2, '[%]')

    # colorbar for diff
    cax3 = h_1_ax(fig, axes[3],ratio=cb_ratio, y_off=cb_y_off_b)
    cb3 = plt.colorbar(diff_pout['mesh'], cax=cax3,
            orientation='horizontal', extend=diff_extend)
    if cb_diff_ticks is not None:
        cb3.set_ticks(cb_diff_ticks)
    right_center_label(cax3, diff_units)
#
#------------------------------------------------------------------------------
#
def plot_scatter_comparison_2(var1, var2, lat_e, lon_e, 
        figsize=(8,6), vmin=None, vmax=None,
        diff_vmin=None, diff_vmax=None,
        rel_pct_vmin=None, rel_pct_vmax=None,
        cb_diff_ticks=None,
        cb_rel_pct_ticks=None,
        extend='neither',
        diff_extend='neither',
        rel_pct_extend='neither',
        left=0.08, right=0.80,
        top=0.95, bottom=0.1,
        hspace=None, wspace=None,
        cl_color='k',
        lw=None,
        region_limit=None,
        xticks=np.arange(-180, 180.0, 60),
        yticks=np.arange(-90, 90.1, 30),
        title_dict={},
        cmap=deepcopy(WhGrYlRd_map),
        diff_cmap=deepcopy(gbcwpry_map),
        rel_pct_cmap=deepcopy(gbcwpry_map),
        units='',
        diff_units='',
        cb_ratio=0.06,
        cb_y_off_t=-0.04,
        cb_y_off_b=-0.07,
        mask_ocean=False,
        ):
    """ Plot var1, var2 and var2 - var1.
    (Yi Wang, 10/21/2020)

    Parameters
    ----------
    var1 : 2-D numpy array
        The first variable.
    var2 : 2-D numpy array
        The second variable.
    lat_e : 2-D numpy array
        Latitude edges.
    lon_e : 2-D numpy array
        Longitude edge.
    figsize : 2-element tuple
        Figure size
    vmin : float or None
        Minimum to be plotted for var1 and var2.
        If vmin is None, it will be set as the valid minimum of var1 
        or var2, which is smaller.
    vmax : float or None
        Similar to vmin, but of maximum
    diff_vmin : float or None
        Minimum to be plotted for var2 - var1.
    diff_vmax : float or None
        Maximum to be plotted for var2 = var1.
    rel_pct_vmin : float or None
        Minimum to be plotted for (var2 - var1) / var1 * 100
    rel_pct_vmax : float or None
        Maximum to be plotted for (var2 - var1) / var1 * 100

    Returns
    -------

    """

    # define figure
    nrow = 2
    ncol = 2
    projPos = list(range(nrow * ncol))
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos,
            countries=True, states=True,
            coastlines=True, mask_ocean=mask_ocean)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']
 
    # get min and max for var1 and var2
    if vmin is None:
        vmin = min(np.nanmin(var1), np.nanmin(var2))
    if vmax is None:
        vmax = max(np.nanmax(var1), np.nanmax(var2))

    # plot var1 and var2
    var_list = [var1, var2]
    varname_list = ['var1', 'var2']
    for i in range(len(var_list)):

        ax = axes[i]
        var = var_list[i]
        title = title_dict.get(varname_list[i], None)
        pout = cartopy_plot_scatter(lon_e, lat_e, c=var, ax=ax, title=title,
                cbar=False, vmin=vmin, vmax=vmax, cmap=deepcopy(cmap))

    # plot (var2 - var1) / var1 * 100
    rel_pct = (var2 - var1) / var1 * 100
    if (rel_pct_vmin is None) and (rel_pct_vmax is None):
        rel_pct_vmax = np.nanmax(np.absolute(rel_pct))
        rel_pct_vmin = -rel_pct_vmax
    ax = axes[2]
    title = title_dict.get('relative_diff', None)
    rel_pct_pout = cartopy_plot_scatter(lon_e, lat_e, c=rel_pct, 
            ax=ax, title=title,
            cbar=False, vmin=rel_pct_vmin, vmax=rel_pct_vmax,
            cmap=deepcopy(rel_pct_cmap))

    # plot var2 - var1
    diff = var2 - var1
    if (diff_vmin is None) and (diff_vmax is None):
        diff_vmax = np.nanmax(np.absolute(diff))
        diff_vmin = -diff_vmax
    ax = axes[3]
    title = title_dict.get('diff', None)
    diff_pout = cartopy_plot_scatter(lon_e, lat_e, c=diff, ax=ax, title=title,
            cbar=False, vmin=diff_vmin, vmax=diff_vmax, 
            cmap=deepcopy(diff_cmap))

    # ticks
    panel_tick_label(axes, ncol, xticks=xticks, yticks=yticks)

    # set limit
    for ax in axes:
        if region_limit is not None:
            ax.set_xlim((region_limit[1],region_limit[3]))
            ax.set_ylim((region_limit[0],region_limit[2]))

    # colorbar for var1 and var2
    cax1 = h_2_ax(fig, axes[0], axes[1],ratio=cb_ratio, y_off=cb_y_off_t)
    cb1 = plt.colorbar(pout['paths'], cax=cax1,
            orientation='horizontal', extend=extend)
    right_center_label(cax1, units)

    # colorbar for (var2 - var1) / var1 * 100
    cax2 = h_1_ax(fig, axes[2],ratio=cb_ratio, y_off=cb_y_off_b)
    cb2 = plt.colorbar(rel_pct_pout['paths'], cax=cax2,
            orientation='horizontal', extend=rel_pct_extend)
    if cb_rel_pct_ticks is not None:
        cb2.set_ticks(cb_rel_pct_ticks)
    right_center_label(cax2, '[%]')

    # colorbar for diff
    cax3 = h_1_ax(fig, axes[3],ratio=cb_ratio, y_off=cb_y_off_b)
    cb3 = plt.colorbar(diff_pout['paths'], cax=cax3,
            orientation='horizontal', extend=diff_extend)
    if cb_diff_ticks is not None:
        cb3.set_ticks(cb_diff_ticks)
    right_center_label(cax3, diff_units)
#
#------------------------------------------------------------------------------
#
