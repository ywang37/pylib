import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PathCollection
#from mpl_toolkits.basemap import Basemap
from copy import deepcopy

def mask_exterior_contiguous_US(ax, m, filename, drawbounds=True, region=[20.0, -130.0, 60.0, -60.0], \
        color='grey'):
    """
    ------------------------------------------------------------------------
    Mask regions that are outside contigious US
    ax: axes
    m: Basemap instance
    filename: US shapefile (cb_2017_us_nation_5m)
    drawbounds: Draw boundaries of shape
    region: [min_lat, min_lon, max_lat, max_lon] of exterior rectangle
    color: fill color
    ------------------------------------------------------------------------
    (Yi Wang, 08/15/2018)
    """

    # read shape file
    m.readshapefile(filename, 'US')

    # m.US[84] is the shapefile vertices of contiguous US
    shape = m.US[84]
    
    Path = mpath.Path

    # interior contiguous US path
    interior_data = []
    interior_data.append( (Path.MOVETO, shape[0]) )
    for i in range(1, len(shape)-1):
        interior_data.append( (Path.LINETO, shape[i]) )
    interior_data.append( (Path.CLOSEPOLY, shape[-1]) )
    in_codes, in_verts = zip(*interior_data)
    interior = mpath.Path(in_verts, in_codes)

    # exterior rectangle
    min_lat = region[0]
    max_lat = region[2]
    min_lon = region[1]
    max_lon = region[3]
    exterior_data = [ \
            (Path.MOVETO,    [min_lon, max_lat]), \
            (Path.LINETO,    [max_lon, max_lat]), \
            (Path.LINETO,    [max_lon, min_lat]), \
            (Path.LINETO,    [min_lon, min_lat]), \
            (Path.CLOSEPOLY, [min_lon, max_lat])  \
            ]
    ex_codes, ex_verts = zip(*exterior_data)
    exterior = mpath.Path(ex_verts, ex_codes)
    exterior.vertices = exterior.vertices[::-1]

    # mask path
    mask_path = mpath.Path(vertices=np.concatenate([exterior.vertices, interior.vertices]), \
            codes=np.concatenate([exterior.codes, interior.codes]), )

    # plot
    collection = PathCollection([mask_path], facecolor=color, zorder=10)
    ax.add_collection(collection)


# test example
if '__main__' == __name__:

    filename = '../data/cb_2017_us_nation_5m/cb_2017_us_nation_5m'

    fig = plt.figure()
    ax = fig.add_subplot(111)

    min_lat = 25.0
    max_lat = 50.0
    min_lon = -125.0
    max_lon = -65.0
    m = Basemap(projection='cyl', resolution='l', llcrnrlat=min_lat, urcrnrlat=max_lat, \
                llcrnrlon=min_lon, urcrnrlon=max_lon, suppress_ticks=True, anchor='C')
    m.drawmeridians(np.arange(-180, 150.1, 10), labels=[0,0,1,0], linewidth=0.0)
    m.drawparallels(np.arange(10,  60.1, 5), labels=[1,0,0,0], linewidth=0.0)

    mask_exterior_contiguous_US(ax, m, filename)

    plt.savefig('mask_exterior_contiguous_US.png', format='png', dpi=300)
    plt.show()





