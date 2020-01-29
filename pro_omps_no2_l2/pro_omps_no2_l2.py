"""
Created on January 28, 2019

@author: Yi Wang
"""

import numpy as np


def QC_OMPS_NO2_L2(NO2, sza, vza, CF, QF_pixel, \
        min_val=-1.0, max_val=10.0, \
        sza_T=75.0, vza_T=65.0, CF_T=0.3):
    """ Quality control for NASA OMPS L2 NO2 (ColumnAmountNO2tropo) 
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
        CloudFraction (Effective cloud fraction)
    QF_pixel : 
        PixelQualityFlag
    min_val : float (default -1.0)
        Threshold for *NO2*, unit is DU
    max_val : float (default 10.0)
        Threshold for *NO2*, unit is DU
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
    val_flag = np.logical_and( NO2 > min_val, NO2 < max_val )

    # geometry
    geo_flag = np.logical_and( sza < sza_T, vza < vza_T )

    # cloud flag
    cloud_flag = ( CF <= CF_T )

    # PixelQualityFlag
    pixel_flag = ( QF_pixel == 0 )

    # final flag
    flag = ( val_flag & geo_flag & cloud_flag & pixel_flag)
    
    return flag
