"""
Created on Feburary 4, 2020

@author: Yi Wang
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#
#------------------------------------------------------------------------------
#
def h_1_ax(fig, ax, x_off=0.0, y_off=-0.1, ratio=0.06, 
        width_add=0.0, height=None):
    """ Add a colobar axes that is under one axes
    for *fig*.

    Parameters
    ----------
    fig : plt.figure() instance
    ax : axes
    x_off : float (default is -0.1)
       Offset from *ax* left
    y_off : float (default is -0.1)
       Offset from *ax* bottom
    ratio : float (default is 0.06)
       The ratio of height to width of the colorbar axes when
       the *height* is not specified
    width_add : float (defaut is 0.0)
       Default width colorbar axes equals that of *ax*.
       *width_add* is used to change the width.
    height : float (default is None)
       The height of the colorbar axes, it override *ratio* when
       *height* is not None, or *height* is determined by *ratio*

    Returns
     -------
    cax : axes
        A colobar axes that is under one axes for *fig*

    """

    # axes poaistion
    pos = ax.get_position()

    # cax
    x0 = pos.x0 + x_off
    y0 = pos.y0 + y_off
    width = pos.width + width_add
    if height is None:
        height = width * ratio
    cax = fig.add_axes([x0, y0, width, height])

    return cax
#
#------------------------------------------------------------------------------
#
def v_1_ax(fig, ax, x_off=0.1, y_off=0.0, ratio=0.06, 
        height_add=0.0, width=None):
    """ Add a colobar axes that is on the right of one axes
    for *fig*.

    Parameters
    ----------
    fig : plt.figure() instance
    ax : axes
    x_off : float (default is 0.1)
       Offset from *ax* left
    y_off : float (default is 0.0)
       Offset from *ax* bottom
    ratio : float (default is 0.06)
       The ratio of width to height of the colorbar axes when
       the *width* is not specified
    width_add : float (defaut is 0.0)
       Default width colorbar axes equals that of *ax*.
       *width_add* is used to change the width.
    width : float (default is None)
       The width of the colorbar axes, it override *ratio* when
       *width* is not None, or *width* is determined by *ratio*

    Returns
     -------
    cax : axes
        A colobar axes that is under one axes for *fig*

    """

    # axes poaistion
    pos = ax.get_position()

    # cax
    x0 = pos.x1 + x_off
    y0 = pos.y0 + y_off
    height = pos.height + height_add
    if width is None:
        width = height * ratio
    cax = fig.add_axes([x0, y0, width, height])

    return cax
#
#------------------------------------------------------------------------------
#
def h_2_ax(fig, ax1, ax2, y_off=-0.1, ratio=0.06, height=None):
    """ Add a colobar axes that is centerd under two axes
    for *fig*.

    Parameters
    ----------
    fig : plt.figure() instance
    ax1 : axes
       left axes
    ax2 : axes
       right axes
    y_off : float (default is -0.1)
       Offset from ax1(ax2) bottom
    ratio : float (default is 0.06)
       The ratio of height to width of the colorbar axes when
       the *height* is not specified
    height : float (default is None)
       The height of the colorbar axes, it override *ratio* when
       *height* is not None, or *height* is determined by *ratio*

    Returns
    -------
    cax : axes
        A colobar axes that is centerd under two axes
        for *fig*

    """

    # axes poaistion
    pos1 = ax1.get_position()
    pos2 = ax2.get_position()

    # cax
    x0 = pos1.x0 + pos1.width / 2.0
    y0 = pos1.y0 + y_off
    width = pos2.x0+pos2.width / 2.0 - (pos1.x0 + pos1.width / 2.0)
    if height is None:
        height = width * ratio

    cax = fig.add_axes([x0, y0, width, height])

    return cax
#
#------------------------------------------------------------------------------
#
def v_2_ax(fig, ax1, ax2, x_off=0.1, y_off=0.0, ratio=0.06,
                height_add=0.0, width=None):
    """ Add a colobar axes that is on the right of two axes
    for *fig*.

    Parameters
    ----------
    fig : plt.figure() instance
    ax1 : axes
        Top axes
    ax2 : axes
        Bottom axes
    x_off : float (default is 0.1)
        Offset from ax1(ax2) right
    y_off : float (default is 0.0)
        Offset from ax2 center
    ratio : float (default is 0.06)
        The ratio of width to height of the colorbar axes when
        the *width* is not specified
    width_add : float (defaut is 0.0)
        Default width colorbar axes equals that of *ax*.
        *width_add* is used to change the width.
    width : float (default is None)
        The width of the colorbar axes, it override *ratio* when
        *width* is not None, or *width* is determined by *ratio*

    Returns
    -------
    cax : axes
        A colobar axes that is under one axes for *fig*

    """

    # axes position
    pos1 = ax1.get_position()
    pos2 = ax2.get_position()

    # cax
    x0 = pos2.x1 + x_off
    y0 = pos2.y0 + pos2.height / 2.0 + y_off
    height = pos1.y0 + pos1.height / 2.0 - \
            (pos2.y0 + pos2.height / 2.0) + height_add
    if width is None:
        width = height * ratio
    cax = fig.add_axes([x0, y0, width, height])

    return cax
#
#------------------------------------------------------------------------------
#
def panel_tick_label(ax_list, ncol, xticks=[], yticks=[],
        xlabel='', ylabel='',
        zero_direction_label=False,
        dateline_direction_label=False, number_format='g',
        degree_symbol=u'\u00B0',
        ):
    """ Add ticks and labels to the right and bottom sides
    of panel plots.
    (Yi Wang, 02/27/2020)

    Parameters
    ----------
    ax_list : list
        All axes are stored in the list.
    ncol : int
        column number of the figure
    xticks : list-like
        xticks
    yticks : list-like
        yticks
    xlabel : str
        xlabel
    ylabel : str
        ylabel
    zero_direction_label :
        Direction label at 0 degree longitude
    dateline_direction_label :
        Direction label at 180 degree longitude
    number_format :
    degree_symbol :

    Returns
    ------

    """

    # xflag
    xflag = [False] * len(ax_list)
    for i in range(ncol):
        xflag[-1-i] = True

    # yflag
    yflag = []
    for i in range(len(ax_list)):
        if (i % ncol) == 0:
            yflag.append(True)
        else:
            yflag.append(False)

    if ax_list is None:
        return xflag, yflag


    # Tick labels
    tick_proj = ['PlateCarree', 'Mercator']

    # Tick format for projection
    lon_formatter = LongitudeFormatter(
            zero_direction_label=zero_direction_label,
            dateline_direction_label=dateline_direction_label,
            number_format=number_format, degree_symbol=degree_symbol)
    lat_formatter = LatitudeFormatter(number_format=number_format,
            degree_symbol=degree_symbol)

    # add ticks and labels
    for i in range(len(ax_list)):

        ax = ax_list[i]

        # projection
        try:
            proj = str(type(ax.projection)).split('.')[-1][0:-2]
        except:
            proj = None

        if xflag[i]:

            # xticks
            if proj in tick_proj:
                ax.set_xticks(xticks, crs=ccrs.PlateCarree())
                ax.xaxis.set_major_formatter(lon_formatter)
            else:
                ax.set_xticks(xticks)

            # xlabel
            ax.set_xlabel(xlabel)

        if yflag[i]:

            # yticks
            if proj in tick_proj:
                ax.set_yticks(yticks, crs=ccrs.PlateCarree())
                ax.yaxis.set_major_formatter(lat_formatter)
            else:
                ax.set_yticks(yticks)

            # ylabel
            ax.set_ylabel(ylabel)

    return xflag, yflag
#
#------------------------------------------------------------------------------
#
def right_center_label(ax, label):
    """ Add *label* to the right center of *ax*.
    This function is usually used to add unit label
    on a colorbar

    Parameters
    ----------
    ax : axes
    label : str
    
    """

    ax.yaxis.set_label_position('right')
    ax.set_ylabel(label, rotation=0, ha='left', va='center')

#
#------------------------------------------------------------------------------
#
def multiFigure(nRow, nCol, ocean_color='grey', **kwargs):
    """
    This function creates multi-figure layout...
	
    Parameter:
    nRow	: row number of the figure
    nCol	: column number of the figure
	
    Optional Paras:
    proj	: projection method (current only allow one projection method
    nUnit	: number of grid per figure
    figsize : figure size
    nGap	: number of grid between each figure
	
    Output:
	
    figure  : figure instance 
    axes	: list of axes instance
    proj	: projection method (current only allow one projection method)
	
    """
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
	
    proj = kwargs.get('proj', ccrs.PlateCarree())
    nUint = kwargs.get('nUint', 1)
    nGap = kwargs.get('nGap', 0)	
    figsize = kwargs.get('figsize', (8, 6))	
    projPos = kwargs.get('projPos', [])
	
    axes = []
    fig = plt.figure(figsize = figsize)

    nRow_grid = nRow * nUint
    nCol_grid = nCol * nUint
    left   = kwargs.get('left', None)
    right  = kwargs.get('right', None)
    top    = kwargs.get('top', None)
    bottom = kwargs.get('bottom', None)
    wspace = kwargs.get('wspace', None)
    hspace = kwargs.get('hspace', None)
    countries = kwargs.get('countries', False)
    states    = kwargs.get('states', False)
    mask_ocean = kwargs.get('mask_ocean', False)
    coastlines = kwargs.get('coastlines', False)
    cl_res = kwargs.get('cl_res', '110m')
    gs = fig.add_gridspec(nRow_grid, nCol_grid, left=left, right=right,
            top=top, bottom=bottom, wspace=wspace, hspace=hspace)

    if states:
        states_provinces = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='50m',
                facecolor='none')

	
    for i in range(nRow):
        for j in range(nCol):
            numFigure = i * nCol + j
            if numFigure in projPos:
                instant = fig.add_subplot(gs[i*nUint:(i + 1)*nUint - nGap,  \
                        j*nUint:(j + 1)*nUint - nGap], \
                        projection=proj)

                # countries
                if countries:
                    instant.add_feature(cfeature.BORDERS)

                # states
                if states:
                    instant.add_feature(states_provinces, 
                            edgecolor='k', linewidth=0.5)

                # mask ocean
                if mask_ocean:
                    instant.add_feature(cfeature.OCEAN, 
                            color=ocean_color, zorder=100)

                # coastlines
                if coastlines:
                    instant.coastlines(zorder=300, resolution=cl_res)

            else:
                instant = fig.add_subplot(gs[i*nUint:(i + 1)*nUint - nGap,  \
                        j*nUint:(j + 1)*nUint - nGap]  )
                
            axes.append(instant)

    out_dict = {}
    out_dict['fig']  = fig
    out_dict['axes'] = axes
    out_dict['proj'] = proj

    return out_dict
#
#------------------------------------------------------------------------------
#
