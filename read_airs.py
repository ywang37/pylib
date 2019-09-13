"""
Created on September 9, 2019

@author: Yi Wang
"""

import copy
from netCDF4 import Dataset

# standard pressure level [hPa]
StdPressureLev = [
        1000. ,  925. ,  850. ,  700. ,  600. ,  500. ,  400. ,
         300. ,  250. ,  200. ,  150. ,  100. ,   70. ,   50. ,
          30. ,   20. ,   15. ,   10. ,    7. ,    5. ,    3. ,
           2. ,    1.5,    1. ,
        ]

def read_AIRS3STD(filename, varnames=[], verbose=False):
    """ read AIRS standard L3 daily product.

    Parameters
    ----------
    filename : str
        AIRS standard L3 daily product file.
    varnames : list
        A list of variable names
    verbose : logical
        Whether or not output more informations.

    Returns
    -------
    out_data : dict
        A dictionary of all varilables.

    """

    if verbose:
        print(' - read_AIRS3STD: reading ' + filename)

    # all variable names
    all_varnames = [
            'Latitude', 'Longitude', 
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

    # add StdPressureLev
    out_data['StdPressureLev'] = copy.deepcopy(StdPressureLev)

    return out_data
