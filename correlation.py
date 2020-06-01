"""
Created on May 27, 2019

@author: Yi Wang
"""

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
        p_vmin=None, p_vmax=None, p_cmap=None,
        countries=True, states=True):
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

    # plot correlation coefficients
    r_plot = cartopy_plot(lon_e, lat_e, r_val, cbar_prop=cbar_prop,
            countries=countries, states=states, cmap=r_cmap,
            vmin=r_vmin, vmax=r_vmax)
    r_plot['cb'].set_label('Linear correlation coefficient')

    # save correlation coefficients plot
    fig_r = fig_dir + name + '_r_' + varname + '.png'
    plt.savefig(fig_r, format='png', dpi=300)

    # plot p value
    p_plot = cartopy_plot(lon_e, lat_e, p_val, cbar_prop=cbar_prop,
            countries=countries, states=states, cmap=p_cmap,
            vmin=p_vmin, vmax=p_vmax)
    p_plot['cb'].set_label('p-value')

    # save p-value
    fig_p = fig_dir + name + '_p_' + varname + '.png'
    plt.savefig(fig_p, format='png', dpi=300)

    # plot correlation coefficients with p-value less than p_thre
    r_val_signi = deepcopy(r_val)
    flag = (p_val < p_thre)
    r_val_signi[np.logical_not(flag)] = np.nan
    r_signi_plot = cartopy_plot(lon_e, lat_e, r_val_signi, 
            cbar_prop=cbar_prop,
            countries=countries, states=states, cmap=r_cmap,
            vmin=r_vmin, vmax=r_vmax)
    r_signi_plot['cb'].set_label('Linear correlation coefficient' + \
            ' (p<{:})'.format(p_thre))

    # save correlation coefficients plot
    fig_r_signi = fig_dir + name + '_r_p' + str(p_thre) + '_' + \
            varname + '.png'
    plt.savefig(fig_r_signi, format='png', dpi=300)




#
#------------------------------------------------------------------------------
#
