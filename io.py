"""
Created on November 6, 2019

@author: Yi Wang
"""

import h5py
from netCDF4 import Dataset
import numpy as np

#
#------------------------------------------------------------------------------
#
def read_hdf5(filename, varnames, verbose=False,
        squeeze=False):
    """ Read HDF5 file
    (ywang, 04/11/20)

    Parameters
    ----------
    filename : str
        HDF5 filename.
    varnames : list
        A list of variable names.
    verbose : logical
        Whether or not output more informations.
    squeeze : logical
        Remove single-dimensional entries from the shape
        of an array.

    Returns
    -------
    out_data : dict
        A dictionary of all variables.
    """

    if verbose:
        print(' - read_hdf5: reading ' + filename)

    # open file
    infile = h5py.File(filename, 'r')

    # read variables
    out_data = {}
    for varn in varnames:
        out_data[varn] = infile[varn][:]

    # close file
    infile.close()

    # Remove single-dimensional entries from the shape
    # of an array
    if squeeze:
        for varn in out_data:
            out_data[varn] = np.squeeze(out_data[varn])

    return out_data
#
#------------------------------------------------------------------------------
#
def read_nc(filename, varnames, verbose=False,
        squeeze=False):
    """ Read netCDF file

    Parameters
    ----------
    filename : str
        netCDF filename.
    varnames : list
        A list of variable names.
    verbose : logical
        Whether or not output more informations.
    squeeze : logical
        Remove single-dimensional entries from the shape
        of an array.

    Returns
    -------
    out_data : dict
        A dictionary of all variables.
    """

    if verbose:
        print(' - read_nc: reading ' + filename)

    # open dataset
    fid = Dataset(filename, 'r')

    # read variables
    out_data = {}
    for varn in varnames:
        tmp = varn.split('/')
        if (len(tmp) == 1):
            out_data[varn] = fid.variables[varn][:]
        else:
            group = fid.groups[tmp[0]]
            for i in range(1, len(tmp)-1):
                group = group.groups[tmp[i]]
            out_data[varn] = group.variables[tmp[-1]][:]

    # close dataset
    fid.close()

    # Remove single-dimensional entries from the shape
    # of an array
    if squeeze:
        for varn in out_data:
            out_data[varn] = np.squeeze(out_data[varn])

    return out_data
#
#------------------------------------------------------------------------------
#
def write_nc(filename, data_dict, units_dict=None, verbose=True):
    """ Write 2-D fields to netCDF file
    (Yi Wang, 02/17/2020)

    Parameters
    ----------
    filename : str
        netCDF filename.
    data_dict : dict
        Variables dictionary
    units_dict : dict
        Unit dictionary
    verbose : logical
        Whether or not output more informations.

    Returns
    -------
    None

    """

    coord_name_list = [
            'Latitude',
            'Latitude_e',
            'Longitude',
            'Longitude_e'
            ]

    if verbose:
        print(' - save_ave: output ' + filename)

    # open file
    nc_f = Dataset(filename, 'w')

    # grid, _e means edge
    Latitude    = data_dict['Latitude']
    Longitude   = data_dict['Longitude']
    Latitude_e  = data_dict['Latitude_e']
    Longitude_e = data_dict['Longitude_e']

    # Dimensions of a netCDF file
    dim_lat = nc_f.createDimension('Latitude',  Latitude.shape[0])
    dim_lon = nc_f.createDimension('Longitude', Latitude.shape[1])
    dim_lat_e = nc_f.createDimension('Latitude_e',  Latitude_e.shape[0])
    dim_lon_e = nc_f.createDimension('Longitude_e', Latitude_e.shape[1])

    # create variables in a netCDF file

    # lat and lon
    Latitude_v = nc_f.createVariable('Latitude', 'f4', 
            ('Latitude', 'Longitude'))
    Longitude_v = nc_f.createVariable('Longitude', 'f4',
            ('Latitude', 'Longitude'))
    Latitude_e_v = nc_f.createVariable('Latitude_e', 'f4',
            ('Latitude_e', 'Longitude_e'))
    Longitude_e_v = nc_f.createVariable('Longitude_e', 'f4',
            ('Latitude_e', 'Longitude_e'))

    # variables
    nc_var_dict = {}
    for varname in data_dict:
        if not (varname in coord_name_list):
            nc_var = nc_f.createVariable(varname, 'f4',
                    ('Latitude', 'Longitude'))
            nc_var_dict[varname] = nc_var

    # write variables

    # lat and lon
    Latitude_v[:]    = Latitude
    Longitude_v[:]   = Longitude
    Latitude_e_v[:]  = Latitude_e
    Longitude_e_v[:] = Longitude_e

    for varname in nc_var_dict:
        nc_var_dict[varname][:] = data_dict[varname]

    # add units
    if units_dict is not None:
        for varname in units_dict:
            nc_var_dict[varname].units = units_dict[varname]

    # close file
    nc_f.close()
#
#------------------------------------------------------------------------------
#
