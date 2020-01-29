"""
Created on Janunary 27, 2019

@author: Yi Wang
"""

import h5py

def read_OMPS_NO2_L2(filename, varnames=[], replace=False, verbose=False):
    """ read OMPS L2 NO2 product

    Parameters
    ----------
    filename : str
        OMPS L2 NO2 product file.
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
        print(' - read_OMPS_NO2_L2: reading ' + filename)

    # all variables names
    all_varnames = [
            '/GeolocationData/Latitude',
            '/GeolocationData/Longitude',
            '/GeolocationData/SolarZenithAngle',
            '/GeolocationData/ViewingZenithAngle',
            '/GeolocationData/SolarAzimuthAngle',
            '/GeolocationData/ViewingAzimuthAngle',
            '/GeolocationData/ImageMidpoint_TAI93',
            '/GeolocationData/GroundPixelQualityFlags',
            '/GeolocationData/InstrumentQualityFlags',
            '/ScienceData/ColumnAmountNO2tropo',
            '/ScienceData/CloudFraction',
            '/ScienceData/PixelQualityFlags',
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

        # get dataset
        var = infile[varname]

        varn = varname.split('/')[-1] # key in out_data

        out_data[varn] = var[:]

    # close dataset
    infile.close()

    return out_data
