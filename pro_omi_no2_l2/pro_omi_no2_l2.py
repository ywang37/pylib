"""
Created on January 14, 2019

@author: Yi Wang
"""

import numpy as np


def QC_OMI_NO2_L2(NO2, sza, vza, CF, XTtrackQ, vcdQ, TerrRef, \
        min_val=-1e30, sza_T=75.0, vza_T=65.0, CF_T=0.3, TerrRef_T=0.3):
    """ Quality control for NASA OMI L2 NO2 (ColumnAmountNO2Trop) 
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
    XTtrackQ : 
        XTrackQualityFlags
    vcdQ :
        VcdQualityFlags
    TerrRef :
        TerrainReflectivity
    min_val : float (default -1e30)
        Threshold for *NO2*
    sza_T : float (default 75.0)
        Threshold for *sza*
    vza_T : float (default 65.0)
        Threshold for *vza*
    CF_T : float (default 0.3)
        Threshold for *CF_T*
    TerrRef_T : float (default 0.3)
        Threshold for *TerrRef*

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

    # XT flag
    XT_flag = np.logical_or(XTtrackQ == 0, XTtrackQ == 255)

    # vcdQ flag
    vcdQ_flag = ( (vcdQ // 2) == 0 )
    for i in range(vcdQ.shape[0]):
        for j in range(vcdQ.shape[1]):
            if ( np.bitwise_and(vcdQ[i,j],2**4) != 0 ):
                vcdQ_flag[i,j] = False

    # TR flag
    TR_flag = ( TerrRef <= TerrRef_T )

    # final flag
    flag = ( val_flag & geo_flag & cloud_flag & XT_flag & vcdQ_flag & TR_flag)
    
    return flag
