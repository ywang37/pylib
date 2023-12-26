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
def write_nc(filename, data_dict, units_dict=None,
        longname_dict=None,
        data_1D_time_dict=None,
        data_1D_time_type_dict={},
        data_3D_time_dict=None, time_type='int',
        verbose=True,
        glo_attr={}):
    """ Write 2-D, 3-D fields to netCDF file
    (Yi Wang, 02/17/2020)

    Parameters
    ----------
    filename : str
        netCDF filename.
    data_dict : dict
        Variables dictionary
    units_dict : dict
        Unit dictionary
    longname_dict : dict
        Long name dictionary
    data_1D_time_dict : dict
        1D varibale dictionary (time,)
    data_1D_time_type_dict : dict
        1D varibale dictionary (time,)
    data_3D_time_dict : dict
        3D varibale dictionary (time, Laititude, Longitude)
    time_type : str
        Data type for time
    verbose : logical
        Whether or not output more informations.
    glo_attr : dict
        Global string attributes

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

    # lat and lon dimensions of a netCDF file
    dim_lat = nc_f.createDimension('Latitude',  Latitude.shape[0])
    dim_lon = nc_f.createDimension('Longitude', Latitude.shape[1])
    dim_lat_e = nc_f.createDimension('Latitude_e',  Latitude_e.shape[0])
    dim_lon_e = nc_f.createDimension('Longitude_e', Latitude_e.shape[1])

    # time dimension
    if (data_1D_time_dict is not None) or \
            (data_3D_time_dict is not None):
        dim_time = nc_f.createDimension('time', None)

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

    # (time, ) varibales
    if data_1D_time_dict is not None:
        for varname in data_1D_time_dict:
            if varname == 'time':
                nc_var = nc_f.createVariable('time', time_type, ('time',))
            else:
                curr_type = data_1D_time_type_dict.get(varname, 'f4')
                nc_var = nc_f.createVariable(varname, curr_type, ('time',))
            nc_var_dict[varname] = nc_var

    # (time, Laititude, Longitude) variables
    if data_3D_time_dict is not None:
        for varname in data_3D_time_dict:
            nc_var = nc_f.createVariable(varname, 'f4', 
                    ('time', 'Latitude', 'Longitude'))
            nc_var_dict[varname] = nc_var

    # write variables

    # lat and lon
    Latitude_v[:]    = Latitude
    Longitude_v[:]   = Longitude
    Latitude_e_v[:]  = Latitude_e
    Longitude_e_v[:] = Longitude_e

    for varname in data_dict:
        if not (varname in coord_name_list):
            nc_var_dict[varname][:] = data_dict[varname]

    if data_1D_time_dict is not None:
        for varname in data_1D_time_dict:
            nc_var_dict[varname][:] = data_1D_time_dict[varname]

    if data_3D_time_dict is not None:
        for varname in data_3D_time_dict:
            nc_var_dict[varname][:] = data_3D_time_dict[varname]

    # add units
    if units_dict is not None:
        for varname in units_dict:
            nc_var_dict[varname].units = units_dict[varname]

    # add longname
    if longname_dict is not None:
        for varname in longname_dict:
            nc_var_dict[varname].longname = longname_dict[varname]

    # add global attributes
    for attrname in glo_attr:
        setattr(nc_f, attrname, glo_attr[attrname])

    # close file
    nc_f.close()
#
#------------------------------------------------------------------------------
#
def write_site_nc(filename, data_dict, units_dict=None,
        data_2D_time_dict=None, time_type='int',
        verbose=True):
    """ Write site data to netCDF file.
    (Yi Wang, 10/20/2020)

    Parameters
    ----------
    filename : str
        netCDF filename.
    data_dict : dict
        Variables dictionary
    units_dict : dict
        Unit dictionary
    data_2D_time_dict : dict
        2D varibale dictionary (time, site)
    time_type : str
        Data type for time
    verbose : logical
        Whether or not output more informations.

    Returns
    -------
    None

    """

    coord_name_list = ['Latitude', 'Longitude']

    str_name_list = ['site_code']

    int_name_list = ['site_i', 'site_j']

    if verbose:
        print(' - save_ave: output ' + filename)

    # open file
    nc_f = Dataset(filename, 'w')

    # grid,
    Latitude    = data_dict['Latitude']
    Longitude   = data_dict['Longitude']

    # site number dimension
    dim_site = nc_f.createDimension('site', Latitude.shape[0])

    # time dimension
    if data_2D_time_dict is not None:
        dim_time = nc_f.createDimension('time', None)

    # create variables in a netCDF file

    # lat and lon
    Latitude_v  = nc_f.createVariable('Latitude',  'f4', ('site',))
    Longitude_v = nc_f.createVariable('Longitude', 'f4', ('site',))

    # variables
    nc_var_dict = {}
    for varname in data_dict:

        if not (varname in coord_name_list+str_name_list+int_name_list):
            nc_var = nc_f.createVariable(varname, 'f4', ('site',))
            nc_var_dict[varname] = nc_var
        elif varname in str_name_list:
            nc_var = nc_f.createVariable(varname, str, ('site',))
            nc_var_dict[varname] = nc_var
        elif varname in int_name_list:
            nc_var = nc_f.createVariable(varname, int, ('site',))
            nc_var_dict[varname] = nc_var

    # (time, site) variables
    if data_2D_time_dict is not None:
        for varname in data_2D_time_dict:
            if varname == 'time':
                nc_var = nc_f.createVariable('time', time_type, ('time',))
            else:
                nc_var = nc_f.createVariable(varname, 'f4',
                        ('time', 'site'))
            nc_var_dict[varname] = nc_var

    # write variables

    # lat and lon
    Latitude_v[:]  = Latitude
    Longitude_v[:] = Longitude

    for varname in data_dict:
        if not (varname in coord_name_list):
            nc_var_dict[varname][:] = data_dict[varname]

    if data_2D_time_dict is not None:
        for varname in data_2D_time_dict:
            nc_var_dict[varname][:] = data_2D_time_dict[varname]

    # add units
    if units_dict is not None:
        for varname in units_dict:
            nc_var_dict[varname].units = units_dict[varname]

    # close file
    nc_f.close()
#
#------------------------------------------------------------------------------
#
