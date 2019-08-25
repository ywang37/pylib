import numpy as np

def calculate_pixel_edge(data):
    """
    using center latitude and longitude to calculate edge
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

