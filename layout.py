"""
Created on Feburary 4, 2020

@author: Yi Wang
"""

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
