import matplotlib.pyplot as plt
import numpy as np

from mylib.cartopy_plot import add_geoaxes, pcolormesh
from mylib.pro_modis.read_modis import get_modis_1km_rgb

#filename = './MYD021KM.A2015344.0725.006.2015344193947.hdf'
filename = '/Dedicated/jwang-data/ywang/project_global_CWA/1km/MODIS_data/\
Gulf_of_Mexico/20181028/MYD021KM.A2018301.1825.061.2018302153109.hdf'

# get data
lat, lon, rgb = get_modis_1km_rgb(filename)

# plot
fig = plt.figure()
xtick=np.arange(-180, 180.1, 5)
ytick=np.arange(-90, 90.1, 5)
ax = add_geoaxes(fig, 111, xtick=xtick, ytick=ytick)
title = '.'.join(filename.split('/')[-1].split('.')[0:4])
pcolormesh(ax, lon, lat, rgb, title=title)

figname = './' + title + '.png'
plt.savefig(figname, format='png', dpi=300)
