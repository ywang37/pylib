"""
Created on September 15, 2019

@author: Yi Wang
"""

import copy
import h5py
from netCDF4 import Dataset


def read_OMI_NO2_L3(filename, varnames=[], replace=False, verbose=False):
    """ read OMI L3 NO2 product

    Parameters
    ----------
    filename : str
        OMI L3 NO2 product file.
    varnames : list
        A list of variable names
    replace : logical
        True: replace all_varnames by varnames
        False: extend varnames to all_varnames
    verbose : logical
        Whether or not output more informations.

    Returns
    -------
    out_data : dict
        A dictionary of all varilables.

    """

    if verbose:
        print(' - read_OMI_NO2_L3: reading ' + filename)

    # all variable names
    all_varnames = [
        '/HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/ColumnAmountNO2TropCloudScreened', 
            ]
    if not replace:
        all_varnames.extend(varnames)
        all_varnames = list(set(all_varnames))
    else:
        all_varnames = varnames

    # open dataset
    infile = h5py.File(filename, 'r')

    # read variables
    out_data = dict()
    for varname in all_varnames:
        varn = varname.split('/')[-1]
        out_data[varn] = infile[varname][:]

    # close dataset
    infile.close()

    return out_data
