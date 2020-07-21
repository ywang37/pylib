"""
Created on July 14, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt

#
#------------------------------------------------------------------------------
#
def plot_time_series_twinx(l_data_dict, r_data_dict, l_varnames,
        r_varnames, xx, l_ratio=False, r_ratio=False,
        l_scale=1.0, r_scale=1.0, xlabel='', l_ylabel='', 
        r_ylabel='', xlim=None, l_ylim=None, r_ylim=None,
        l_color_dict={}, r_color_dict={},
        l_ylabel_color='C0', r_ylabel_color='C1',
        l_ytick_color='C0', r_ytick_color='C1',
        l_label_dict={}, r_label_dict={}, legend=True,
        fig=None, ax=None):
    """ Plot time series that share the same x-axis.
    (ywang, 07/14/20)

    Parameters
    ----------
    l_data_dict : dict
        Data for the left y-axis. Keys are elements in the
        *l_varnames*, and values are 1-D numpy array.
    r_data_dict : dict
        Similar to *l_data_dict*, but for the right y-axis.
    l_varnames : list
        Elements are keys in the *l_data_dict*
    r_varnames : list
        Elements are keys in the *r_data_dict*
    xx : 1-D numpy array
        Data for the x-axis.
    l_ratio : bool
        If True, each time series in the *l_data_dict* is
        divided by their first elements in the time series.
    r_ratio : bool
        Similar to *l_ratio*, but for *r_data_dict*
    l_scale : float
        Data in *l_data_dict* (maybe normalized (l_ratio)) are 
        multiplied by *l_scale*
    r_scale : float
        Data in *r_data_dict* (maybe normalized (r_ratio)) are
        multiplied by *r_scale*
    xlabel : str
        x-axis label.
    l_ylabel : str
        Left y-axis label.
    r_ylabel : str
        Right y-axis label.
    xlim : None or list like [minv, maxv]
        x-axis limit.
    l_ylim : None or list like [minv, maxv]
        Left y-axis limit.
    r_ylim : None or list like [minv, maxv]
        Right y-axis limit.
    l_color_dict : dict
        Colors for *l_data_dict*
    r_color_dict : dict
        Color for *r_data_dict*
    l_ylabel_color : str
        Left y-axis label color
    r_ylabel_color : str
        Right y-axis label color
    l_ytick_color : str
        Left y-axis tick color
    r_ytick_color : str
        Right y-axis tick color
    l_label_dict : dict
        Labels for *l_data_dict*
    r_label_dict : dict
        Labels for *r_data_dict*
    legend : bool
        Plot legend or not.
    fig : None or plt.figure() instance.
        If both *fig* and *ax* are None, a plt.figure() instance
        will be created.
    ax : None or axes instance.
        If None, an axes instance will be created.

    Returns
    -------
    out_dict : dict
        Keys are 'fig', 'ax', and 'ax_t'

    """

    # axes
    if ax is None:
        if fig is None:
            fig = plt.figure()
        ax = fig.add_subplot(111)

    out_dict = {}
    out_dict['fig'] = fig
    out_dict['ax']  = ax

    ax.set_xlabel(xlabel)

    if xlim is not None:
        ax.set_xlim(xlim)

    #-------------
    # Left y-axis
    #-------------

    l_line_list = []
    for l_varn in l_varnames:

        l_yy = deepcopy(l_data_dict[l_varn])
        if l_ratio:
            l_yy = l_yy / l_yy[0]

        l_label = l_label_dict.get(l_varn, '')
        l_color = l_color_dict.get(l_varn, 'C0')

        l_line, = ax.plot(xx, l_yy * l_scale, marker='o', label=l_label,
                c=l_color)
        l_line_list.append(l_line)

    ax.set_ylabel(l_ylabel)

    if l_ylim is not None:
        ax.set_ylim(l_ylim)

    ax.yaxis.label.set_color(l_ylabel_color)

    ax.tick_params(axis='y', colors=l_ytick_color)

    #--------------
    # Right y-axis
    #--------------

    ax_t = ax.twinx()
    out_dict['ax_t'] = ax_t

    r_line_list = []
    for r_varn in r_varnames:

        r_yy = deepcopy(r_data_dict[r_varn])
        if r_ratio:
            r_yy = r_yy / r_yy[0]

        r_label = r_label_dict.get(r_varn, '')
        r_color = r_color_dict.get(r_varn, 'C1')

        r_line, = ax_t.plot(xx, r_yy * r_scale, marker='o', label=r_label,
                c=r_color)
        r_line_list.append(r_line)

    ax_t.set_ylabel(r_ylabel)
    if r_ylim is not None:
        ax_t.set_ylim(r_ylim)

    ax_t.yaxis.label.set_color(r_ylabel_color)

    ax_t.tick_params(axis='y', colors=r_ytick_color)

    #
    plt.subplots_adjust(left=0.20, right=0.80, bottom=0.25, top=0.75)

    if legend:
        line_list = l_line_list + r_line_list
        ax.legend(line_list, [l.get_label() for l in line_list],
                loc='best')

    return out_dict
#
#------------------------------------------------------------------------------
#
