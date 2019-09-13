"""
Created on September 9, 2019

@author: Yi Wang
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset

from mylib.cartopy_plot import add_geoaxes, pcolormesh
from mylib.pro_satellite import calculate_pixel_edge2
from mylib.read_airs import read_AIRS3STD

def plot_AIRS3STD_lev(in_data, varnames, levels,
        scale_dict={}, vmin_dict={}, vmax_dict={},
        cbar_label_dict={}, title_dict={},
        fig_dir='./', fig_pre=None, fig_post=None,
        verbose=True):
    """ plot AIRS standard L3 daily product at specific
    levels. 

    Parameters
    ----------
    in_data : dict or str
        dict includes all variables for plot.
        str is filename
    varnames : list
        A list of variable names
    levels : list
        pressure levels
    scale_dict : dict
        Some variables to be scaled in plot
    vmin_dict : dict
    vmax_dict : dict
    cbar_label_dict : dict
    title_dict : dict
    fig_dir : str
        directory to save figures
    fig_pre : str or None
        prefix of figure name
    fig_post : str or None
        postfix of figure name
    verbose : logical
        Whether or not output more informations.

    Returns
    -------
    out_data : dict
        A dictionary of all varilables.

    """

    valid_min = 0.0

    # read data if input is str
    if isinstance(in_data, str):
        in_data = read_AIRS3STD(in_data, varnames, verbose=verbose)

    # get standard pressure level [hPa]
    StdPressureLev = in_data['StdPressureLev']

    # loop variable
    for varn in varnames:

        # get variable
        varv = in_data[varn]

        # loop pressure
        for lev in levels:

            ind = StdPressureLev.index(lev)
            pres = StdPressureLev[ind]

            if verbose:
                print(' - plot_AIRS3STD_lev: plot '
                    + varn 
                    + ' at level {}, pressure is {} hPa.'.format(ind, pres)
                    )
            
            # plot figure
            fig = plt.figure(figsize=(8,5))
            ax = add_geoaxes(fig, 111)
            lat = in_data['Latitude']
            lon = in_data['Longitude']
            lat_e, lon_e = calculate_pixel_edge2(lat, lon)
            varv_lev = varv[ind,:,:]
            sl = scale_dict.get(varn, 1.0)
            varv_lev = varv_lev * sl
            vmin = vmin_dict.get(varn, None)
            vmax = vmax_dict.get(varn, None)
            title = fig_post + ' ' + str(int(pres)) + 'hPa ' \
                    + title_dict.get(varn, varn)
            mesh = pcolormesh(ax, lon_e, lat_e, varv_lev, 
                    valid_min=valid_min, title=title, 
                    vmin=vmin, vmax=vmax)
            plt.subplots_adjust(bottom=0.2)

            # add colorbar
            cax = fig.add_axes([0.2, 0.13, 0.6, 0.03])
            cbar = plt.colorbar(mesh, cax=cax, orientation='horizontal')
            cbar_label = cbar_label_dict.get(varn, '')
            cbar.set_label(cbar_label)

            # save figure
            fig_name = varn + '_' + str(int(pres)) + 'hPa'
            if fig_pre is not None:
                fig_name = fig_pre + '_' + fig_name
            if fig_post is not None:
                fig_name = fig_name + '_' + fig_post
            fig_name = fig_dir + fig_name + '.png'
            plt.savefig(fig_name, format='png', dpi=300)
            plt.close()

def plot_AIRS3STD_col(in_data, varnames,
        scale_dict={}, vmin_dict={}, vmax_dict={},
        cbar_label_dict={}, title_dict={},
        fig_dir='./', fig_pre=None, fig_post=None,
        verbose=True):
    """ plot AIRS standard L3 daily vertical column density.

    Parameters
    ----------
    in_data : dict or str
        dict includes all variables for plot.
        str is filename
    varnames : list
        A list of variable names
    scale_dict : dict
        Some variables to be scaled in plot
    vmin_dict : dict
    vmax_dict : dict
    cbar_label_dict : dict
    title_dict : dict
    fig_dir : str
        directory to save figures
    fig_pre : str or None
        prefix of figure name
    fig_post : str or None
        postfix of figure name
    verbose : logical
        Whether or not output more informations.

    Returns
    -------
    out_data : dict
        A dictionary of all varilables.

    """

    valid_min = 0.0

    # read data if input is str
    if isinstance(in_data, str):
        in_data = read_AIRS3STD(in_data, varnames, verbose=verbose)

    # get standard pressure level [hPa]
    StdPressureLev = in_data['StdPressureLev']

    # loop variable
    for varn in varnames:

        if verbose:
            print(' - plot_AIRS3STD_lev: plot '
                + varn + ' vertical column density'
                )
            
        # plot figure
        fig = plt.figure(figsize=(8,5))
        ax = add_geoaxes(fig, 111)
        lat = in_data['Latitude']
        lon = in_data['Longitude']
        lat_e, lon_e = calculate_pixel_edge2(lat, lon)
        varv = in_data[varn]
        sl = scale_dict.get(varn, 1.0)
        varv_sl = varv * sl
        vmin = vmin_dict.get(varn, None)
        vmax = vmax_dict.get(varn, None)
        title = fig_post + ' ' + title_dict.get(varn, varn)
        mesh = pcolormesh(ax, lon_e, lat_e, varv_sl, 
                valid_min=valid_min, title=title, 
                vmin=vmin, vmax=vmax)
        plt.subplots_adjust(bottom=0.2)

        # add colorbar
        cax = fig.add_axes([0.2, 0.13, 0.6, 0.03])
        cbar = plt.colorbar(mesh, cax=cax, orientation='horizontal')
        cbar_label = cbar_label_dict.get(varn, '')
        cbar.set_label(cbar_label)

        # save figure
        fig_name = varn 
        if fig_pre is not None:
            fig_name = fig_pre + '_' + fig_name
        if fig_post is not None:
            fig_name = fig_name + '_' + fig_post
        fig_name = fig_dir + fig_name + '.png'
        plt.savefig(fig_name, format='png', dpi=300)
        plt.close()




















