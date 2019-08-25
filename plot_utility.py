import matplotlib.pyplot as plt
import numpy as np


def plot_box(m, region, **kwargs):
    """ Plot a rectangle on map

    Parameters
    ----------
    m: Basemap instance

    region : tuple-like
    Define the region to be ploted.
    (min_lat, min_lon, max_lat, max)

    """

    num = 50

    # left side
    x_l = np.array([region[1]]*num)
    y_l = np.array(np.linspace(region[0], region[2], num=num))
    m.plot(x_l, y_l, latlon=True, **kwargs)

    # right side
    x_r = np.array([region[3]]*num)
    y_r = np.array(np.linspace(region[0], region[2], num=num))
    m.plot(x_r, y_r, latlon=True, **kwargs)

    # top side
    x_t = np.array(np.linspace(region[1], region[3], num=num))
    y_t = np.array([region[2]]*num)
    m.plot(x_t, y_t, latlon=True, **kwargs)

    # bottom side
    x_b = np.array(np.linspace(region[1], region[3], num=num))
    y_b = np.array([region[0]]*num)
    m.plot(x_b, y_b, latlon=True, **kwargs)


