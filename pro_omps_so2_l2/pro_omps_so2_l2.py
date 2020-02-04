"""
Created on January 28, 2019

@author: Yi Wang
"""

import numpy as np


def QC_OMPS_SO2_L2(SO2, sza, vza, RCF, QF_pixel, \
        min_val=-5.0, max_val=15.0, \
        sza_T=75.0, vza_T=65.0, RCF_T=0.2):
    """ Quality control for NASA OMPS L2 SO2 (ColumnAmountSO2_PBL) 
    product.

    Parameters
    ----------
    SO2 : numpy array
        ColumnAmountSO2_PBL
    sza : numpy array
        SolarZenithAngle
    vza : numpy array
        ViewZenithAngle
    RCF : numpy array
        RadiativeCloudFraction
    QF_pixel : 
        PixelQualityFlag
    min_val : float (default -5.0)
        Threshold for *SO2*, unit is DU
    max_val : float (default 15.0)
        Threshold for *SO2*, unit is DU
    sza_T : float (default 75.0)
        Threshold for *sza*
    vza_T : float (default 65.0)
        Threshold for *vza*
    RCF_T : float (default 0.3)
        Threshold for *RCF_T*

    Returns
    -------
    flag : logical array
        True: valid
        False: invalid
    
    """

    # value
    val_flag = np.logical_and( SO2 > min_val, SO2 < max_val )

    # geometry
    geo_flag = np.logical_and( sza < sza_T, vza < vza_T )

    # cloud flag
    cloud_flag = ( RCF <= RCF_T )

    # PixelQualityFlag
    pixel_flag = ( QF_pixel == 0 )

    # final flag
    flag = ( val_flag & geo_flag & cloud_flag & pixel_flag)
    
    return flag
