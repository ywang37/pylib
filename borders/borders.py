"""
Created on September 2, 2020

@author: Yi Wang
"""

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt

#
#------------------------------------------------------------------------------
#
def border_US(ax, res='110m', plot_stat=True):
    """ Plot US border.
    (ywang, 06/01/2021)

    Parameters
    ----------
    ax : GeoAxes
        axes
    res : str
        Resolution. '10m', '50m', or '110m'
    plot_stat : bool
        Plot states or not

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
            ax.add_geometries(country.geometry, ccrs.PlateCarree(),
                    facecolor='none', edgecolor='k', lw=1)

    # read states
    if plot_stat:
        shpfstat = shpreader.natural_earth(resolution=res,
                category='cultural',
                name='admin_1_states_provinces_shp')
                # 'admin_1_states_provinces_shp' has country boundries
                # 'admin_1_states_provinces_lines' does not have
                # country boundries, and use 'adm0_a3'
        readerstat = shpreader.Reader(shpfstat)
        states = readerstat.records()

    # plot states
    for state in states:
        if state.attributes['sr_adm0_a3'] == 'USA':
            ax.add_geometries(state.geometry, ccrs.PlateCarree(),
                    facecolor='none', edgecolor='k', 
                    lw=0.5)

    return None
#
#------------------------------------------------------------------------------
#

if '__main__'== __name__:

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    border_US(ax)
    ax.set_extent([-180, -65, 0, 90])

    plt.show()
