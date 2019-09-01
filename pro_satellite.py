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


