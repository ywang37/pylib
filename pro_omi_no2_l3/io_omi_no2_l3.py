"""
Created on September 16, 2019

@author: Yi Wang
"""

from netCDF4 import Dataset

from mylib.io import read_hdf5

#
#------------------------------------------------------------------------------
#
def read_OMI_NO2_L3(filename, verbose=True):
    """ Read OMI L3 NO2
    (ywang, 04/11/20)

    Parameters
    ----------
    filename : str
        OMI L3 NO2 data file.
    verbose : logical
        Whether or not output more informations.

    Returns
    -------
    out_data : dict
        A dictionary of all variables.
    """

    print(' - read_OMI_NO2_L3: reading ' + filename)

    # variables
    short_varnames = [ \
            'ColumnAmountNO2TropCloudScreened',
            ]
    data_path = '/HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/'
    varnames = []
    for varn in short_varnames:
        varnames.append(data_path + varn)

    # read data
    out_dict = read_hdf5(filename, varnames)

    # change to short variable names
    for varn in varnames:
        tmp = varn.split('/')[-1]
        out_dict[tmp] = out_dict.pop(varn)

    return out_dict
#
#------------------------------------------------------------------------------
#
def output_month_OMI_NO2_L3(out_dir, month, NO2, verbose=False):
    """ Output monthly mean of OMI L3 NO2

    Parameters
    ----------
    out_dir : str
        Output directory
    month : str
        Month, used for generate filename.
    NO2 :
        Monthly mean of OMI L3 NO2
    """ 

    filename = out_dir + 'OMI_NO2_' + month + '_monthly.nc'

    if verbose:
        print(' - output_month_OMI_NO2_L3: writing' + filename)

    # open a file
    nc_f = Dataset(filename, 'w', format='NETCDF4')

    # Dimensions of a netCDF file
    dim = NO2.shape
    nlat_c = dim[0]
    nlon_c = dim[1]

    dim_lat_c = nc_f.createDimension('lat_c', nlat_c)
    dim_lon_c = nc_f.createDimension('lon_c', nlon_c)

    # Create variables in a netCDF file
    NO2_v = nc_f.createVariable('NO2', 'f4', ('lat_c', 'lon_c'))

    # longname
    NO2_v.longname = 'Monthly mean ColumnAmountNO2CloudScreened'

    # units
    NO2_v.units = 'molec/cm2'

    # write variable
    NO2_v[:,:] = NO2

    # close a file
    nc_f.close()

def output_multi_year_month_OMI_NO2_L3(out_dir, start_year, end_year, 
        month, NO2, verbose=False):
    """ Output multi-year mean of monthly  OMI L3 NO2.

    Parameters
    ----------
    out_dir : str
        Output directory
    start_year : str
        Used to generate filename and longname
    end_year : str
        Used to generate filename and longname
    month : str
        Used to generate filename and longname
    NO2 : 2D numpy array
        multi-year mean of monthly  OMI L3 NO2
    verbose : logical
        output more information

    """

    filename = out_dir + 'OMI_NO2_' + start_year + '-' + end_year \
            + '_' + month + '_monthly.nc'

    if verbose:
        print(' - output_multi_year_month_OMI_NO2_L3: output ' + filename)
        
    # open a file
    nc_f = Dataset(filename, 'w', format='NETCDF4')

    # Dimensions of a netCDF file
    dim = NO2.shape
    nlat_c = dim[0]
    nlon_c = dim[1]

    dim_lat_c = nc_f.createDimension('lat_c', nlat_c)
    dim_lon_c = nc_f.createDimension('lon_c', nlon_c)

    # Create variables in a netCDF file
    NO2_v = nc_f.createVariable('NO2', 'f4', ('lat_c', 'lon_c'))

    # longname
    NO2_v.longname = month + ' monthly mean ColumnAmountNO2CloudScreened (' \
            + start_year + '-' + end_year + ')'

    # units
    NO2_v.units = 'molec/cm2'

    # write variable
    NO2_v[:,:] = NO2

    # close a file
    nc_f.close()

def read_month_OMI_NO2_L3(filename, verbose=False):
    """ Read monthly mean of OMI L3 NO2

    Parameters
    ----------
    filename : str
        monthly mean of OMI L3 NO2 file
    verbose : logical
        output more information

    Returns
    -------
    NO2 : 2D numpy array
        Monthly mean of OMI L3 NO2

    """

    if verbose:
        print(' - read_month_OMI_NO2_L3: read ' + filename)

    nc_f = Dataset(filename)

    NO2 = nc_f.variables['NO2'][:]

    nc_f.close()

    return NO2
    
def read_multi_year_month_OMI_NO2_L3(in_dir, start_year, end_year, 
        month, verbose=False):
    """ Output multi-year mean of monthly  OMI L3 NO2.

    Parameters
    ----------
    out_dir : str
        Output directory
    start_year : str
        Used to generate filename and longname
    end_year : str
        Used to generate filename and longname
    month : str
        Used to generate filename and longname
    verbose : logical
        output more information

    Returns
    -------
    NO2 : 2D numpy array
        multi-year mean of monthly  OMI L3 NO2

    """

    filename = in_dir + 'OMI_NO2_' + start_year + '-' + end_year \
            + '_' + month + '_monthly.nc'

    if verbose:
        print(' - read_multi_year_month_OMI_NO2_L3: read ' + filename)

    nc_f = Dataset(filename)

    NO2 = nc_f.variables['NO2'][:]

    nc_f.close()

    return NO2








