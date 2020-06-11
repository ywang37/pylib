"""
Created on November 5, 2019

@author: Yi Wang
"""

from copy import deepcopy
import numpy as np

#
#------------------------------------------------------------------------------
#
def tropomi_no2_qc(NO2, min_val=None, max_val=None,
        CRFNW=None, CRFNW_thre=0.3, 
        QA=None, QA_thre=0.5, SZA=None, SZA_thre=70.0,
        VZA=None, VZA_thre=70.0):
    """ Quality control for TROPOMI
    nitrogendioxide_tropospheric_column.
    (ywang, 06/01/2020)

    Parameters
    ----------
    NO2 : numpy array 
        nitrogendioxide_tropospheric_column
    min_val : None or 
        Minimum NO2
    max_val : None or
        Maximum NO2
    CRFNW : numpy array or None
        cloud_radiance_fraction_nitrogendioxide_window
        (Cloud radiance fraction at 440 nm for NO2 retrieval)
    CRFNW_thre : float 
       CRFNW threshold
    QA: numpy array or None
        qa_value
        (A continuous quality descriptor, varying between 
        0 (no data) and 1 (full quality data). Recommend to 
        ignore data with qa_value < 0.5)
    QA_thre : float o
        QA threshod
    SZA : numpy array or None
        Solar zenith angle
    SZA_thre : float
        Solar zenith angle threshold
    VZA : numpy array or None
        View zenith angle
    VZA_thre : float
        View zenith angle threshold

    Returns
    -------
    out_dict : dictornary
        NO2_QC : numpy array
            Data after quality control. Data with low quality
            are replace by np.nan
        nan_flag : numpy logical array

    """

    # 
    nan_flag = np.full_like(NO2, False)

    # min_val
    if min_val is not None:
        nan_flag = np.logical_or(nan_flag, NO2 < min_val)

    # max_val
    if max_val is not None:
        nan_flag = np.logical_or(nan_flag, NO2 > max_val)

    # CRFNW
    if CRFNW is not None:
        nan_flag = np.logical_or(nan_flag, CRFNW > CRFNW_thre)

    # QA
    if QA is not None:
        nan_flag = np.logical_or(nan_flag, QA < QA_thre)

    # SZA
    if SZA is not None:
        nan_flag = np.logical_or(nan_flag, SZA > SZA_thre)

    # VZA 
    if VZA is not None:
        nan_flag = np.logical_or(nan_flag, VZA > VZA_thre)

    # NO2_QC
    NO2_QC = deepcopy(NO2)
    NO2_QC[nan_flag] = np.nan

    # out_dict
    out_dict = {}
    out_dict['NO2_QC'] = NO2_QC
    out_dict['nan_flag'] = nan_flag

    return out_dict
#
#------------------------------------------------------------------------------
#
