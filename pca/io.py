"""
Created on April 12, 2020

@author: Yi Wang
"""

from netCDF4 import Dataset
import numpy as np

#
#------------------------------------------------------------------------------
#
def write_1D_PCA_nc(filename, data_dict, units_dict=None, verbose=True):
    """ Write 1-D PCA result to netCDF file
    (Yi Wang, 03/17/2022)

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
            'Longitude',
            ]

    if verbose:
        print(' - save_ave: output ' + filename)

    # open file
    nc_f = Dataset(filename, 'w')

    # latitude and longitude
    Latitude    = data_dict['Latitude']
    Longitude   = data_dict['Longitude']

    # other data
    ori_data   = data_dict['ori_data']
    ave        = data_dict['ave']
    components = data_dict['components']
    coeffs     = data_dict['coeffs']
    explained_variance_ratio = data_dict['explained_variance_ratio']
    explained_variance = data_dict['explained_variance']

    # Dimensions of a netCDF file
    dim_samp = nc_f.createDimension('n_samples',    ori_data.shape[0]  )
    dim_site = nc_f.createDimension('n_sites',      ori_data.shape[1]  )
    dim_comp = nc_f.createDimension('n_components', components.shape[0])

    # create variables in a netCDF file

    # lat and lon
    Latitude_v = nc_f.createVariable('Latitude',   'f4', ('n_sites',))
    Longitude_v = nc_f.createVariable('Longitude', 'f4', ('n_sites',))

    # variables
    ori_data_v = nc_f.createVariable('ori_data', 'f4', 
            ('n_samples', 'n_sites'))
    ave_v = nc_f.createVariable('ave', 'f4', ('n_sites',))
    components_v = nc_f.createVariable('components', 'f4',
            ('n_sites', 'n_components'))
    coeffs_v = nc_f.createVariable('coeffs', 'f4',
            ('n_samples', 'n_components'))
    explained_variance_ratio_v = \
            nc_f.createVariable('explained_variance_ratio', 'f4',
            ('n_components',))
    explained_variance_v = nc_f.createVariable('explained_variance', 'f4',
            ('n_components',))

    # write variables

    # lat and lon
    Latitude_v[:]    = Latitude
    Longitude_v[:]   = Longitude

    # data
    ori_data_v[:]                 = ori_data
    ave_v[:]                      = ave
    components_v[:]               = components
    coeffs_v[:]                   = coeffs
    explained_variance_ratio_v[:] = explained_variance_ratio
    explained_variance_v[:]       = explained_variance

    # add units
    if units_dict is not None:
        for varname in units_dict:
            nc_var_dict[varname].units = units_dict[varname]

    # close file
    nc_f.close()
#
#------------------------------------------------------------------------------
#
def write_2D_PCA_nc(filename, data_dict, units_dict=None, verbose=True):
    """ Write 2-D PCA result to netCDF file
    (Yi Wang, 04/12/2020)

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

    # other data
    ori_data   = data_dict['ori_data']
    ave        = data_dict['ave']
    components = data_dict['components']
    coeffs     = data_dict['coeffs']
    explained_variance_ratio = data_dict['explained_variance_ratio']
    explained_variance = data_dict['explained_variance']

    # Dimensions of a netCDF file
    dim_lat = nc_f.createDimension('Latitude',  Latitude.shape[0])
    dim_lon = nc_f.createDimension('Longitude', Latitude.shape[1])
    dim_lat_e = nc_f.createDimension('Latitude_e',  Latitude_e.shape[0])
    dim_lon_e = nc_f.createDimension('Longitude_e', Latitude_e.shape[1])
    dim_comp = nc_f.createDimension('n_components', components.shape[2])
    dim_samp = nc_f.createDimension('n_samples', ori_data.shape[2])

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
    ori_data_v = nc_f.createVariable('ori_data', 'f4',
            ('Latitude', 'Longitude', 'n_samples'))
    ave_v = nc_f.createVariable('ave', 'f4',
            ('Latitude', 'Longitude'))
    components_v = nc_f.createVariable('components', 'f4',
            ('Latitude', 'Longitude', 'n_components'))
    coeffs_v = nc_f.createVariable('coeffs', 'f4',
            ('n_samples', 'n_components'))
    explained_variance_ratio_v = \
            nc_f.createVariable('explained_variance_ratio', 'f4',
            ('n_components',))
    explained_variance_v = nc_f.createVariable('explained_variance', 'f4',
            ('n_components',))

    # write variables

    # lat and lon
    Latitude_v[:]    = Latitude
    Longitude_v[:]   = Longitude
    Latitude_e_v[:]  = Latitude_e
    Longitude_e_v[:] = Longitude_e

    # data
    ori_data_v[:]                 = ori_data
    ave_v[:]                      = ave
    components_v[:]               = components
    coeffs_v[:]                   = coeffs
    explained_variance_ratio_v[:] = explained_variance_ratio
    explained_variance_v[:]       = explained_variance

    # add units
    if units_dict is not None:
        for varname in units_dict:
            nc_var_dict[varname].units = units_dict[varname]

    # close file
    nc_f.close()
#
#------------------------------------------------------------------------------
#
