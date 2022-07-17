"""
Created on March 11, 2019

@author: Yi Wang
"""

import numpy as np

def QC_GEMS_NO2_L2(NO2, sza, vza, CF,\
        min_val=-1e30, sza_T=75.0, vza_T=65.0, CF_T=0.3):
    """ Quality control for GEMS L2 NO2 (ColumnAmountNO2Trop)
    product.

    Parameters
    ----------
    NO2 : numpy array
        ColumnAmountNO2Trop
    sza : numpy array
        SolarZenithAngle
    vza : numpy array
        ViewZenithAngle
    CF : numpy array
        CloudFraction
    min_val : float (default -1e30)
        Threshold for *NO2*
    sza_T : float (default 75.0)
        Threshold for *sza*
    vza_T : float (default 65.0)
        Threshold for *vza*
    CF_T : float (default 0.3)
        Threshold for *CF_T*

    Returns
    -------
    flag : logical array
        True: valid
        False: invalid

    """

    # value
    val_flag = ( NO2 > min_val )

    # geometry
    geo_flag = np.logical_and( sza < sza_T, vza < vza_T )

    # cloud flag
    cloud_flag = ( CF <= CF_T )

    # final flag
    flag = ( val_flag & geo_flag & cloud_flag )

    return flag
