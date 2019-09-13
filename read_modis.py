"""
Created on September 11, 2019

@author: Yi Wang
"""

from geotiepoints import modis5kmto1km
from netCDF4 import Dataset
import numpy as np

def get_modis_1km_rgb(filename, verbose=False):
    """ get RGB reflectance from MxD021KM file.

    Parameters
    ----------
    filename : str
        MxD021KM file.

    Returns
    -------
    rgb : 3D numpy array
        R, G, and B channel TOA reflectance
        (without solar zenith angel correction.)

    """

    varname_list = [\
            'EV_250_Aggr1km_RefSB', \
            'EV_500_Aggr1km_RefSB', \
            'EV_500_Aggr1km_RefSB', \
            ]

    ind_list = [0, 1, 0]

    if verbose:
        print(' - get_modis_1km_rgb: reading ' + filename)

    # open dataset
    fid = Dataset(filename, 'r')

    # get RGB
    rgb = []
    for i in range(len(varname_list)):

        varname = varname_list[i]
        ind = ind_list[i]

        var = fid.variables[varname]
        var.set_auto_maskandscale(False)

        scale_factor = var.reflectance_scales[ind]
        offset       = var.reflectance_offsets[ind]

        var = var
        var = var[ind,:,:]
        var = (var - offset) * scale_factor

        var[var<0.0] = 0.0
        var[var>1.0] = 1.0

        if i == 0:
            rgb = np.zeros((var.shape[0],var.shape[1],3))

        rgb[:,:,i] = var

    # get latitude and lonitude
    lat = fid.variables['Latitude'][:]
    lon = fid.variables['Longitude'][:]
    lon, lat = modis5kmto1km(lon, lat)


    # close dataset
    fid.close()

    return lat, lon, rgb

def read_mxd04(filename, varnames=[], verbose=False):
    """ Read MxD04_L2 product.

    Parameters
    ----------
    filename : str
        MxD04_L2 file.
    varnames : list
        A list of variable names.
    verbose : logical
        Whether or not output more informations.

    Returns
    -------
    out_data : dict
        A dictionary of all variables.
    """

    if verbose:
        print(' - read_mxd04: reading ' + filename)

    # all variable names
    all_varnames = [
            'Latitude', 'Longitude',
            'Effective_Optical_Depth_Average_Ocean',
            'Mean_Reflectance_Ocean',
            'Optical_Depth_Land_And_Ocean',
            ]
    all_varnames.extend(varnames)
    all_varnames = list(set(all_varnames))

    # open dataset
    fid = Dataset(filename, 'r')

    # read variables
    out_data = {}
    for varname in all_varnames:
        out_data[varname] = fid.variables[varname][:]

    # close dataset
    fid.close()

    return out_data







