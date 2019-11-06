"""
Created on November 5, 2019

@author: Yi Wang
"""

#
#------------------------------------------------------------------------------
#

def tropomi_no2_qc(NO2, CRFNW=None, CRFNW_thre=None, 
        QA=None, QA_thre=None):
    """ Quality control for TROPOMI
    nitrogendioxide_tropospheric_column.

    Parameters
    ----------
    NO2 : numpy array 
        nitrogendioxide_tropospheric_column
    CRFNW : numpy array or None
        cloud_radiance_fraction_nitrogendioxide_window
        (Cloud radiance fraction at 440 nm for NO2 retrieval)
    CRFNW_thre : float or None
       CRFNW threshold
    QA: numpy array
        qa_value
        (A continuous quality descriptor, varying between 
        0 (no data) and 1 (full quality data). Recommend to 
        ignore data with qa_value < 0.5)
    QA_thre : float or None
        QA threshod

    Returns
    -------
    NO2_mask : numpy masked array
        Data after quality control

    """

#
#------------------------------------------------------------------------------
#

