"""
Created on March 10, 2022

@author: Yi Wang
"""

#from netCDF4 import Dataset
from mylib.io import read_nc

def read_GEMS_NO2_L2(filename, varnames=[], replace=False, verbose=False):
    """ read GEMS L2 NO2 product

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
        print(' - read_GEMS_NO2_L2: reading ' + filename)

    # all variables names
    all_varnames = [
            'Geolocation Fields/Latitude',
            'Geolocation Fields/Longitude',
            'Geolocation Fields/SolarZenithAngle',
            'Geolocation Fields/ViewingZenithAngle',
            'Data Fields/CloudFraction',
            'Data Fields/ColumnAmountNO2Trop',
            ]
    if not replace:
        all_varnames.extend(varnames)
        all_varnames = list(set(all_varnames))
    else:
        all_varnames = varnames

    # read data
    out_data = read_nc(filename, all_varnames, verbose=False, squeeze=False)

    return out_data
