"""
Created on May 27, 2019

@author: Yi Wang
"""

import cartopy.feature as cfeature
import cartopy.crs as ccrs
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

from mylib.cartopy_plot import cartopy_plot
from mylib.io import read_nc

#
#------------------------------------------------------------------------------
#
def plot_pearsonr_map(filename, fig_dir, varname, name='',
        p_thre=0.05, r_vmin=None, r_vmax=None,
        r_cmap=plt.get_cmap('seismic'),
        r_signi_cmap=plt.get_cmap('seismic'),
        r_signi_vmin=None, r_signi_vmax=None,
        r_signi_min_valid=None,
        p_vmin=None, p_vmax=None, p_cmap=None,
        signi_ocean_color=None,
        countries=True, states=True,
        cl_res='110m',
        cl_color=None,
        lw=None,
        equator=False, NH_label='', SH_label=''):
    """ Plot pearson correlation coefficent.
    (ywang, 05/27/20)

    Parameters
    ----------
    filename : str
        netCDF file of trend_analysis results.
    fig_dir : str
        Directory to save figures.
    varname : str
        variable name
    name : str
        Prefix of figure names.
    p_thre : float
        p-value threshod.
    signi_ocean_color : None or str
        If None, the *signi_ocean_color*  color to
        mask r_signi map if it is not None

    """

    # directory
    if fig_dir[-1] != '/':
        fig_dir = fig_dir + '/'

    # variables to be read
    varnames = ['Latitude_e', 'Longitude_e', \
            'r_'+varname, 'p_'+varname]

    # read data
    data = read_nc(filename, varnames, verbose=True)

    # get latitude and longitude edges
    lat_e = data['Latitude_e']
    lon_e = data['Longitude_e']

    # get correlation coefficients and p-value
    r_val = data['r_'+varname]
    p_val = data['p_'+varname]

    cbar_prop = {}
    cbar_prop['orientation'] = 'horizontal'

    # equator line
    eqr_c = 'white'
    eqr_ls = '--'
    eqr_ls_lw = 2

    # plot correlation coefficients
    r_plot = cartopy_plot(lon_e, lat_e, r_val, cbar_prop=cbar_prop,
            countries=countries, states=states, cmap=r_cmap,
            vmin=r_vmin, vmax=r_vmax)
    r_plot['cb'].set_label('Linear correlation coefficient')
    if equator:
        r_plot['ax'].plot([-180, 180], [0, 0], color=eqr_c,
                linestyle=eqr_ls, lw=eqr_ls_lw, transform=ccrs.PlateCarree())
        r_plot['ax'].text(-175, 5, NH_label, color=eqr_c, ha='left',
                va='bottom', transform=ccrs.Geodetic())
        r_plot['ax'].text(-175, -5, SH_label, color=eqr_c, ha='left',
                va='top', transform=ccrs.Geodetic())


    # save correlation coefficients plot
    fig_r = fig_dir + name + '_r_' + varname + '.png'
    plt.savefig(fig_r, format='png', dpi=300)

    # plot p value
    p_plot = cartopy_plot(lon_e, lat_e, p_val, cbar_prop=cbar_prop,
            countries=countries, states=states, cmap=p_cmap,
            vmin=p_vmin, vmax=p_vmax)
    p_plot['cb'].set_label('p-value')
    if equator:
        p_plot['ax'].plot([-180, 180], [0, 0], color=eqr_c,
                linestyle=eqr_ls, lw=eqr_ls_lw, transform=ccrs.PlateCarree())
        p_plot['ax'].text(-175, 5, NH_label, color=eqr_c, ha='left',
                va='bottom', transform=ccrs.Geodetic())
        p_plot['ax'].text(-175, -5, SH_label, color=eqr_c, ha='left',
                va='top', transform=ccrs.Geodetic())

    # save p-value
    fig_p = fig_dir + name + '_p_' + varname + '.png'
    plt.savefig(fig_p, format='png', dpi=300)

    # plot correlation coefficients with p-value less than p_thre
    r_val_signi = deepcopy(r_val)
    flag = (p_val < p_thre)
    r_val_signi[np.logical_not(flag)] = np.nan
    if r_signi_min_valid is not None:
        r_val_signi[r_val_signi < r_signi_min_valid] = np.nan
    r_signi_plot = cartopy_plot(lon_e, lat_e, r_val_signi, 
            cbar_prop=cbar_prop,
            countries=countries, states=states, cmap=r_signi_cmap,
            vmin=r_signi_vmin, vmax=r_signi_vmax)
    r_signi_plot['cb'].set_label('Linear correlation coefficient' + \
            ' (p<{:})'.format(p_thre))
    if equator:
        r_signi_plot['ax'].plot([-180, 180], [0, 0], color=eqr_c,
                linestyle=eqr_ls, lw=eqr_ls_lw, transform=ccrs.PlateCarree(),
                zorder=200)
        r_signi_plot['ax'].text(-175, 5, NH_label, color=eqr_c, ha='left',
                va='bottom', transform=ccrs.Geodetic(), zorder=200)
        r_signi_plot['ax'].text(-175, -5, SH_label, color=eqr_c, ha='left',
                va='top', transform=ccrs.Geodetic(), zorder=200)
    if signi_ocean_color is not None:
        r_signi_plot['ax'].add_feature(cfeature.OCEAN, 
                color=signi_ocean_color, zorder=100)
        r_signi_plot['ax'].coastlines(resolution=cl_res, 
                color=cl_color, lw=lw, zorder=300)

    # save correlation coefficients plot
    fig_r_signi = fig_dir + name + '_r_p' + str(p_thre) + '_' + \
            varname + '.png'
    plt.savefig(fig_r_signi, format='png', dpi=300)




#
#------------------------------------------------------------------------------
#
