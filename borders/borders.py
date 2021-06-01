"""
Created on September 2, 2020

@author: Yi Wang
"""

import cartopy.io.shapereader as shpreader

#
#------------------------------------------------------------------------------
#
def border_US(ax, res='110m'):
    """ Plot US border.
    (ywang, 06/01/2021)

    Parameters
    ----------
    ax : GeoAxes
        axes
    res : str
        Resolution. '10m', '50m', or '110m'

    Returns
    -------
    None

    """

    # read countries
    shpfcoun = shpreader.natural_earth(resolution=res,
            category='cultural',
            name='admin_0_countries')
    readercoun = shpreader.Reader(shpfcoun)
    countries = readercoun.records()

    # plot US
    for country in countries:
        if country.attributes['ADM0_A3'] == 'USA':
            ax.add_geometries()



    return None
#
#------------------------------------------------------------------------------
#
