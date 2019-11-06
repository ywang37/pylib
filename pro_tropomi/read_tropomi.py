"""
Created on November 11, 2019

@author: Yi Wang
"""

from mylib.io import read_nc

#
#------------------------------------------------------------------------------
#

def read_tropomi_no2(filename, varnames=[], replace=False, verbose=False,
        squeeze=True):
    """ Read TROPOMI NO2 product.

    Parameters
    ----------
    filename : str
        S5P_RPRO_L2__NO2*.nc file.
    varnames : list
        A list of variable name.
    replace : logical (default False)
        If True, replcae, else append
    verbose : logical (default False)
        Whether or not output more information.

    Returns
    -------
    out_data : dict
        A dictionary of variables.

    """

    if verbose:
        print(' - read_tropomi_no2: reading ' + filename)

    # all variables names
    all_varnames = [
            'PRODUCT/latitude', 
            'PRODUCT/longutide',
            'PRODUCT/nitrogendioxide_tropospheric_column',
            'PRODUCT/qa_value',
            'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/\
cloud_radiance_fraction_nitrogendioxide_window'
            ]
    if replace:
        all_varnames = varnames
    else:
        all_varnames.extend(varnames)

    # remove duplicated variables
    all_varnames = set(all_varnames)
    all_varnames = list(all_varnames)

    # read data
    out_data = read_nc(filename, all_varnames, squeeze=squeeze)

    return out_data

#
#------------------------------------------------------------------------------
#

