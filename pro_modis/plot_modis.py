"""
Created on Dec 4, 2019

@author: Yi Wang
"""

import datetime
import numpy as np

from mylib.cartopy_plot import cartopy_plot
from mylib.pro_modis.read_modis import get_modis_1km_rgb

def plot_modis_1km_rgb(mxd021km_file, ax=None, fig=None,
        title='date', xtick=None, ytick=None, cl_res=None,
        region_limit=None,
        verbose=False):
    """ Read MxD021KM file to plot RGB image.
    (Yi Wang, 12/04/2019)

    Parameters
    ----------
    mxd021km_file : str
        MxD021KM file.
    ax : GeoAxes or None (default)
        Create a GeoAxes if ax is None.
    fig : plt.figure() or None (default)
        Create a plt.figure() if both fig and ax are None
    title : None or str
        None: no title
        'date': date calcuated from filename
        other str: title

    Returns
    -------

    """

    # ticks
    if xtick is None:
        xtick = np.arange(-180, 180.1, 10)
    if ytick is None:
        ytick = np.arange(-90, 90.1, 5)

    # get RGB, lat, and lon
    lat, lon, rgb = get_modis_1km_rgb(mxd021km_file)

    # title
    if title == 'date':
        tmp = mxd021km_file.split('/')[-1].split('.')
        year = tmp[1][1:5]
        day_of_year = int(tmp[1][5:8])
        curr_date_d = \
                datetime.datetime.strptime(year+'-01-01', '%Y-%m-%d') \
                + datetime.timedelta(days=(day_of_year-1))
        ymd = str(curr_date_d)[0:10]
        title = ymd

    # plot RGB
    cartopy_plot(lon, lat, rgb, ax=ax, fig=fig,
            title=title, xtick=xtick, ytick=ytick,
            cl_res=cl_res, region_limit=region_limit)
