"""
Created on April 12, 2020

@author: Yi Wang
"""

import cartopy.feature as cfeature
from copy import deepcopy
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

from mylib.cartopy_plot import cartopy_plot
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.layout import multiFigure, h_1_ax, h_2_ax, panel_tick_label

#
#------------------------------------------------------------------------------
#
def plot_2D_components_and_coeffs(components, lat_e, lon_e,
        coeffs, explained_variance_ratio,
        fig_dir, name='', n_components=5,
        xticks=np.arange(-180.0, 180.1, 60.0),
        yticks=np.arange(-90.0, 90.1, 30.0),
        region_limit=[-90.0, -180.0, 90.0, 180.0],
        time_ticks=None,
        rotation=30.0,
        reverse = [],
        ):
    """ plot 2D PCA spatial pattern
    (ywang, 04/12/20)

    Parameters
    ----------
    components : shape (n_lat, n_lon, n_components)
        PCA modes
    lat_e : shape (n_lat+1, n_lon+1)
        Latitude edges
    lon_e : shape (n_lat+1, n_lon+1)
        Longitude edges
    coeffs : shape (n_samples, n_components)
        PCA coefficients
    explained_variance_ratio : shape (n_components, )
        Explained variance ratio
    fig_dir : str
        Directory to save figures
    n_components : int
        The first *n_component* modes and coeffs will be plotted.
    time_ticks : tuple
        The first element is tick positions, and the second
        element is tick labels.
    reverse : list
        Mode and time series to be reversed. 1-based

    """

    # directory
    if fig_dir[-1] != '/':
        fig_dir = fig_dir + '/'

    # n_components
    n_components = min([n_components, components.shape[2], coeffs.shape[1]])

    for ic in range(n_components):

        # begin plot 
        nrow = 2
        ncol = 1
        figsize = (8, 6)
        projPos = [0]
        layout_dict = multiFigure(nrow, ncol,
                figsize=figsize, projPos=projPos)
        fig  = layout_dict['fig']
        axes = layout_dict['axes']

        # reverse
        if (ic+1) in reverse:
            sign = -1
        else:
            sign = 1

        # plot mode
        mode = components[:,:,ic] * sign
        vmax = np.max(np.absolute(mode))
        vmin = -vmax
        pout = cartopy_plot(lon_e, lat_e, mode, 
                ax=axes[0], vmin=vmin, vmax=vmax,
                cmap=deepcopy(gbcwpry_map), cbar=True)
        states_provinces = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='50m',
                facecolor='none')
        axes[0].add_feature(cfeature.BORDERS)
        axes[0].add_feature(states_provinces, edgecolor='k', linewidth=0.5)
        axes[0].add_feature(cfeature.COASTLINE, zorder=200)
        axes[0].set_xlim((region_limit[1],region_limit[3]))
        axes[0].set_ylim((region_limit[0],region_limit[2]))

        # add title
        y_p = 1.03
        left_title = 'Mode {}'.format(ic+1)
        axes[0].text(0, y_p, left_title, 
                ha='left', va='bottom', transform=axes[0].transAxes)
        right_title = '{:.1f}%'.format(explained_variance_ratio[ic]*100)
        axes[0].text(1, y_p, right_title, 
                ha='right', va='bottom', transform=axes[0].transAxes)

        # ticks
        panel_tick_label(axes[0:1], ncol, xticks=xticks, yticks=yticks)

        
        # plot coeff
        coeff = coeffs[:,ic] * sign
        x = range(len(coeff))
        axes[1].plot(x, coeff, '-o', markersize=3)
        axes[1].set_xlabel('Month')

        if time_ticks is not None:
            axes[1].set_xticks(time_ticks[0])
            axes[1].set_xticklabels(time_ticks[1], rotation=rotation)

        # save figure
        if (name != ''):
            figname = fig_dir + name + '_mode{}'.format(ic+1) + '.png'
        else:
            figname = fig_dir + 'mode{}'.format(ic+1) + '.png'
        plt.savefig(figname, format='png', dpi=300)

#
#------------------------------------------------------------------------------
#
def plot_explained_variance_ratio(explained_variance_ratio, ax=None):
    """ Plot variance contribution.
    (ywang, 04/12/20)

    Parameters
    ----------
    explained_variance_ratio : 1D array
        Explained variance ratio

    """

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.subplots_adjust(left=0.25, right=0.75, 
                bottom=0.25, top=0.75)
    
    # accumulated contribution
    accum = np.zeros_like(explained_variance_ratio)
    for i in range(len(accum)):
        accum[i] = np.sum(explained_variance_ratio[0:i+1])

    # begin plot
    x = np.array(range(len(explained_variance_ratio))) + 1
    ax.bar(x, explained_variance_ratio*100, label='Each')
    ax.plot(x, accum*100, '-o', c='C1', label='Cumulative')
    ax.set_ylim([0,100])
    ax.set_xlabel('Mode')
    ax.set_ylabel('Variance contributions [%]')
    plt.legend(loc='best')
#
#------------------------------------------------------------------------------
#
def plot_ave(ave, lat_e, lon_e,
        xticks=np.arange(-180.0, 180.1, 60.0),
        yticks=np.arange(-90.0, 90.1, 30.0),
        region_limit=[-90.0, -180.0, 90.0, 180.0],
        vmin=None, vmax=None,
        units='',
        ):
    """ plot average
    (ywang, 04/12/20)

    Parameters
    ----------
    ave : 2D array
        PCA modes
    lat_e : shape (n_lat+1, n_lon+1)
        Latitude edges
    lon_e : shape (n_lat+1, n_lon+1)
        Longitude edges

    """


    # begin plot 
    nrow = 2
    ncol = 1
    figsize = (8, 6)
    projPos = [0]
    layout_dict = multiFigure(nrow, ncol,
            figsize=figsize, projPos=projPos)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']

    # plot average
    pout = cartopy_plot(lon_e, lat_e, ave, 
            ax=axes[0], vmin=vmin, vmax=vmax,
            cmap=deepcopy(WhGrYlRd_map), cbar=True)
    pout['cb'].set_label(units)
    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')
    panel_tick_label(axes[0:1], ncol, xticks=xticks, yticks=yticks)
    axes[0].add_feature(cfeature.BORDERS)
    axes[0].add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    axes[0].add_feature(cfeature.COASTLINE, zorder=200)
    axes[0].set_xlim((region_limit[1],region_limit[3]))
    axes[0].set_ylim((region_limit[0],region_limit[2]))

#
#------------------------------------------------------------------------------
#
