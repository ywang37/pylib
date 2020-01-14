import numpy as np

def calculate_pixel_edge(data):
    """
    using center latitude(longitude) to calculate edge
    """

    dim = data.shape
    num_l = dim[0]
    num_p = dim[1]

    data_e = np.zeros((num_l+1,num_p+1))


    data_e_tmp = (data[0:num_l-1,:] + data[1:num_l,:]) * 0.5
    data_e_tmp = (data_e_tmp[:,0:num_p-1] + data_e_tmp[:,1:num_p]) * 0.5
    data_e[1:num_l,1:num_p] = np.array(data_e_tmp)
    data_e[0,:]  = 2.0 * data_e[1,:]  - data_e[2,:]
    data_e[-1,:] = 2.0 * data_e[-2,:] - data_e[-3,:]
    data_e[:,0]  = 2.0 * data_e[:,1]  - data_e[:,2]
    data_e[:,-1] = 2.0 * data_e[:,-2] - data_e[:,-3]

    return data_e

def calculate_pixel_edge2(lat_c, lon_c):
    """
    Using latitude center and longotude center to calculate
    latitude edge and longitude edge.
    """

    lat_e = calculate_pixel_edge(lat_c)
    lon_e = calculate_pixel_edge(lon_c)

    return lat_e, lon_e

def get_granule_corner(lat, lon):
    """Get the latitude and longitude of granule four corners.
    """

    if (lat.shape != lon.shape):
        print('get_granule_corner: the dimensions of lat and lon are \
inconsistent')

    NY = lat.shape[0]
    NX = lat.shape[1]

    # indexes of the four corners
    corn_index = ( (0, 0), (0, NX-1), (NY-1, NX-1), (NY-1, 0) )

    # latitude and longitude of the four corns
    corn = []
    for cind in range(len(corn_index)):

        i = corn_index[cind][0]
        j = corn_index[cind][1]

        corn_lat_lon = (lat[i,j], lon[i,j])

        corn.append(corn_lat_lon)
    
    return tuple(corn)

def get_granule_center(lat, lon):
    """Get the latitude and longitude of granule four center
    """

    if (lat.shape != lon.shape):
        print('get_granule_corner: the dimensions of lat and lon are \
inconsistent')

    NY = lat.shape[0]
    NX = lat.shape[1]

    # center at Y direction
    if ( (NY % 2) == 0 ):

        i2 = NY // 2
        i1 = i2 - 1

        lat_y = (lat[i1,:] + lat[i2,:]) * 0.5
        lon_y = (lon[i1,:] + lon[i2,:]) * 0.5

    else:

        i = NY // 2

        lat_y = lat[i,:]
        lon_y = lon[i,:]

    # center at X direction
    if ( (NX % 2) == 0 ):

        j2 = NX // 2
        j1 = j2 - 1

        lat_c = (lat_y[j1] + lat_y[j2]) * 0.5
        lon_c = (lon_y[j1] + lon_y[j2]) * 0.5

    else:

        j = NX // 2

        lat_c = lat_y[j]
        lon_c = lon_y[j]

    return (lat_c, lon_c)

def bytescale(data, cmin=None, cmax=None, high=255, low=0):
    """
    Byte scales an array (image).
    Byte scaling means converting the input image to uint8 dtype and scaling
    the range to ``(low, high)`` (default 0-255).
    If the input image already has dtype uint8, no scaling is done.
    Parameters
    ----------
    data : ndarray
        PIL image data array.
    cmin : scalar, optional
        Bias scaling of small values. Default is ``data.min()``.
    cmax : scalar, optional
        Bias scaling of large values. Default is ``data.max()``.
    high : scalar, optional
        Scale max value to `high`.  Default is 255.
    low : scalar, optional
        Scale min value to `low`.  Default is 0.
    Returns
    -------
    img_array : uint8 ndarray
        The byte-scaled array.
    Examples
    --------
    >>> from scipy.misc import bytescale
    >>> img = np.array([[ 91.06794177,   3.39058326,  84.4221549 ],
    ...                 [ 73.88003259,  80.91433048,   4.88878881],
    ...                 [ 51.53875334,  34.45808177,  27.5873488 ]])
    >>> bytescale(img)
    array([[255,   0, 236],
           [205, 225,   4],
           [140,  90,  70]], dtype=uint8)
    >>> bytescale(img, high=200, low=100)
    array([[200, 100, 192],
           [180, 188, 102],
           [155, 135, 128]], dtype=uint8)
    >>> bytescale(img, cmin=0, cmax=255)
    array([[91,  3, 84],
           [74, 81,  5],
           [52, 34, 28]], dtype=uint8)
    """
    if data.dtype == np.uint8:
        return data

    if high < low:
        raise ValueError("`high` should be larger than `low`.")

    if cmin is None:
        cmin = data.min()
    if cmax is None:
        cmax = data.max()

    cscale = cmax - cmin
    if cscale < 0:
        raise ValueError("`cmax` should be larger than `cmin`.")
    elif cscale == 0:
        cscale = 1

    scale = float(high - low) / cscale
    bytedata = (data * 1.0 - cmin) * scale + 0.4999
    bytedata[bytedata > high] = high
    bytedata[bytedata < 0] = 0
    return np.cast[np.uint8](bytedata) + np.cast[np.uint8](low)
    
def scale_image(image, x=None, y=None):
    """ color enhancement
    """

    if x is None:
        x = np.array([0,  30,  60, 120, 190, 255], dtype=np.uint8)

    if y is None:
        y = np.array([0, 110, 160, 210, 240, 255], dtype=np.uint8)

    # bytescale
    image_bs = bytescale(image)

    scaled = np.zeros_like(image_bs, dtype=np.uint8)
    for i in range(len(x)-1):
        x1 = x[i]
        x2 = x[i+1]
        y1 = y[i]
        y2 = y[i+1]
        m = (y2 - y1) / float(x2 - x1)
        b = y2 - (m *x2)
        mask = ((image_bs >= x1) & (image_bs < x2))
        scaled = scaled + mask * np.asarray(m * image_bs + b, dtype=np.uint8)

    mask = image_bs >= x2
    scaled = scaled + (mask * 255)

    scaled = scaled.astype(float)

    scaled = scaled / 255.0

    return scaled

def scale_image_multi_c(image, x=None, y=None):
    """ color enhancement of mutliple channnels.
    """

    scaled = np.zeros_like(image)

    for i in range(image.shape[-1]):
        scaled[:,:,i] = scale_image(image[:,:,i], x=x, y=y)

    return scaled






