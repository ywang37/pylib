"""
Created on Janunary 14, 2019

@author: Yi Wang
"""

import h5py
#import numpy as np
#from netCDF4 import Dataset

def read_OMI_NO2_L2(filename, varnames=[], replace=False, verbose=False):
    """ read OMI L2 NO2 product

    Parameters
    ----------
    filename : str
        OMI L2 NO2 product file.
    varnames : list
        A list of variable names
    replace : logical
        True: replace all_varnames_dict by varnames
        False: update varnames_dict to all_varnames
    verbose : logical
        Whether or not output more informations.

    Returns
    -------
    out_data : dict
        A dictionary of all varilables.

    """

    if verbose:
        print(' - read_OMI_NO2_L2: reading ' + filename)

    # all variables names
    all_varnames = [
            '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Latitude',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Longitude',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/SolarZenithAngle',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/ViewingZenithAngle',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/SolarAzimuthAngle',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/ViewingAzimuthAngle',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Time',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/XTrackQualityFlags',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/VcdQualityFlags',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/AmfTrop',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Trop',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWeight',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWtPressure',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/TerrainReflectivity',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/CloudFraction',
            ]
    if not replace:
        all_varnames.extend(varnames)
        all_varnames = list(set(all_varnames))
    else:
        all_varnames = varnames

    # variables that are retrieved as
    # (value - Offset) * ScaleFactor
    all_scl_off_varnames = [
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/TerrainReflectivity',
            '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/CloudFraction',
            ]

    # open dataset
    infile = h5py.File(filename, 'r')

    # read variables
    out_data = dict()
    for varname in all_varnames:

        # get dataset
        var = infile[varname]

        varn = varname.split('/')[-1] # key in out_data
        if varname in all_scl_off_varnames:
            out_data[varn] =  (var[:] - var.attrs['Offset']) \
                    * var.attrs['ScaleFactor']
        else:
            out_data[varn] = var[:]

    # close dataset
    infile.close()

    return out_data
