"""
Created on January 23, 2020

@author: Yi Wang
"""

import cartopy.crs as ccrs
import numpy as np
import os

def add_China_province(ax, lw=None):
    """ Add China province borderlines

    Paramaters
    ----------
    ax : GeoAxes

    """

    filename = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + \
            '/CN-border-La.dat'
    with open(filename) as src:
        context = src.read()
        blocks = [cnt for cnt in context.split('>') if len(cnt) > 0]
        borders = \
                [np.fromstring(block, dtype=float, sep=' ') \
                for block in blocks]

    if lw is None:
        lw = 1

    # Plot border lines
    for line in borders:
        ax.plot(line[0::2], line[1::2], '-', lw=lw, color='k',
                transform=ccrs.Geodetic())
