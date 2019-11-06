"""
Created on November 6, 2019

@author: Yi Wang
"""

from netCDF4 import Dataset
import numpy as np

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
                group = group.groups[tmp[0]]
            out_data[varn] = group.variables[tmp[-1]][:]

    # close dataset
    fid.close()

    # Remove single-dimensional entries from the shape
    # of an array
    if squeeze:
        for varn in out_data:
            out_data[varn] = np.squeeze(out_data[varn])

    return out_data







