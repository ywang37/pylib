"""
Created on September 17, 2019

@author: Yi Wang
"""

import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

from mylib.cartopy_plot import add_geoaxes, pcolormesh
from mylib.grid_utility import generate_grid, get_center_index
import mylib.pro_omi_no2_l3.io_omi_no2_l3 as io

def plot_NO2_monthly_mean(ax, lat_e, lon_e, NO2):
    """ 
    """

    mesh = pcolormesh(ax, lon_e, lat_e, NO2)


def plot_NO2_monthly_anomaly(ax, lat_e, lon_e, NO2_anomaly):
    """
    """

    mesh = pcolormesh(ax, lon_e, lat_e, NO2_anomaly)


def plot_NO2_monthly_relative_anomaly(ax, lat_e, lon_e, NO2_rel_ano):
    """
    """

    NO2_rel_ano_percent = NO2_rel_ano
    mesh = pcolormesh(ax, lon_e, lat_e, NO2_rel_ano_percent)


def plot_NO2_combine(dir_mean, dir_month,
        start_year, end_year, year, month, verbose=False,
        region_limit=[-90.0, -180.0, 90.0, 180.0], dir_fig=None,
        vmin_mean=0.0, vmax_mean=0.5, vmin_ano=-0.1, vmax_ano=0.1,
        vmin_ano_per=-50.0, vmax_ano_per=50.0):
    """ Plot mutli-year mean of monthly NO2, monthly NO2
    anomaly and monthly NO2 relative anomaly.

    Parameters
    ----------
    dir_mean : str
        Directory of mutli-year mean of monthly NO2
    dir_month : str
        Directory of monthly NO2
    start_year : str
        Used to generate filename and longname
    end_year : str
        Used to generate filename and longname
    year : str
        Used to generate filename and longname
    month : str
        Used to generate filename and longname
    verbose : logical
        output more information
    region_limit : list
        The region for plot. [min_lat, min_lon, max_lat, max_lon]
    dir_fig : str or None
        If it is str, save figure to the directory

    """

    # used to convert unit
    toDU = 2.69e16

    # colorbar parameters
    pad = 0.01
    aspect = 15

    # read monthly data
    file_month = dir_month + 'OMI_NO2_' + year + '-' +  month \
            + '_monthly.nc'
    NO2_month = io.read_month_OMI_NO2_L3(file_month, verbose=verbose)
    NO2_month = NO2_month / toDU

    # read mutli-year mean of monthly NO2
    NO2_mean = io.read_multi_year_month_OMI_NO2_L3(dir_mean, 
            start_year, end_year, month, verbose=verbose)
    NO2_mean = NO2_mean / toDU

    # NO2 anoanmly
    NO2_ano = NO2_month - NO2_mean

    # NO2 realtive anoamly
    NO2_ano_percent = NO2_ano / NO2_mean * 100.0

    # generate grid
    nlat_c, nlon_c = 720, 1440
    lat_e, lon_e, lat_c, lon_c = generate_grid(nlat_c, nlon_c)

    # get region data
    i1 = get_center_index(lat_e, region_limit[0])
    i2 = get_center_index(lat_e, region_limit[2])
    j1 = get_center_index(lon_e, region_limit[1])
    j2 = get_center_index(lon_e, region_limit[3])
    NO2_mean = NO2_mean[i1:i2+1,j1:j2+1]
    NO2_ano = NO2_ano[i1:i2+1,j1:j2+1]
    NO2_ano_percent = NO2_ano_percent[i1:i2+1,j1:j2+1]
    lat_e = lat_e[i1:i2+2]
    lon_e = lon_e[j1:j2+2]

    # begin plot
    #fig = plt.figure(figsize=(4, 8))
    fig = plt.figure()
    xtick=np.arange(-180, 180.1, 20)
    ytick=np.arange(-90, 90.1, 10)
    ax_list = []

    # mutli-year mean of monthly NO2
    title1 = start_year + '-' + end_year + ', ' + month \
            + r', NO$_2$'
    ax1 = add_geoaxes(fig, 311, xtick=[], ytick=ytick, title=title1)
    ax_list.append(ax1)
    mesh1 = pcolormesh(ax1, lon_e, lat_e, NO2_mean,
            vmin=vmin_mean, vmax=vmax_mean)
    cbar1 = fig.colorbar(mesh1, ax=ax1, extend='max',
            pad=pad, aspect=aspect)
    cbar1.set_label(r'[DU]')

    # anomaly
    title2 = year + ', ' + month + ', anomaly'
    ax2 = add_geoaxes(fig, 312, xtick=[], ytick=ytick, title=title2)
    ax_list.append(ax2)
    mesh2 = pcolormesh(ax2, lon_e, lat_e, NO2_ano, 
            cmap=plt.get_cmap('seismic'),
            vmin=vmin_ano, vmax=vmax_ano)
    cbar2 = fig.colorbar(mesh2, ax=ax2, extend='both',
            pad=pad, aspect=aspect)
    cbar2.set_label(r'[DU]')

    # relative anomaly
    title3 = year + ', ' + month + ', relative anomaly'
    ax3 = add_geoaxes(fig, 313, xtick=xtick, ytick=ytick, title=title3)
    ax_list.append(ax3)
    mesh3 = pcolormesh(ax3, lon_e, lat_e, NO2_ano_percent, 
            cmap=plt.get_cmap('seismic'),
            vmin=vmin_ano_per, vmax=vmax_ano_per)
    cbar3 = fig.colorbar(mesh3, ax=ax3, extend='both', 
            pad=pad, aspect=aspect)
    cbar3.set_label(r'[%]')

    # set axes
    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none'
            )
    for ax in ax_list:
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, zorder=200)
        ax.add_feature(cfeature.OCEAN, color='w', zorder=100)
        ax.set_xlim((region_limit[1],region_limit[3]))
        ax.set_ylim((region_limit[0],region_limit[2]))

    plt.subplots_adjust(hspace=0.3)

    # save figure
    if dir_fig is not None:
        figname = dir_fig + 'NO2_anomaly_' + year + '-' + month \
                + '_from_' + start_year + '-' + end_year + '.png'
    plt.savefig(figname, format='png', dpi=300)


    return fig
















