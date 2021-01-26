"""
Created on Sep 4, 2020

@author: Yi Wang
"""

from scipy import stats
import numpy as np

#
#------------------------------------------------------------------------------
#
def inter_1km_to_hkm(var_1km):
    """ Interpolate variables at 1 km resoltion into
    half km resolution
    (Yi Wang, 09/30/2020)

    Parameters
    ----------
    var_1km : 2-D numpy array
        Satellite variable at the resolution of 1 km

    Returns
    -------
    var_hkm : 2-D numpy array
        Satellite variable at the resolution of half km
    var_hkm_edge : 2-D numpy array
         Satellite variable at the resolution of half km at edge
  
    """

    # dimension
    dim_1km = var_1km.shape
    dim_hkm = (dim_1km[0]*2, dim_1km[1]*2)

    # define half km arrays
    var_hkm = np.zeros(dim_hkm)
    var_hkm_edge = np.zeros((dim_hkm[0]+1, dim_hkm[1]+1))

    # 1 km variable is at half km even edge
    var_hkm_edge[1:-1:2, 1:-1:2] = var_1km

    # the rest half km edge (not including the outer edges)
    var_hkm_edge[2:-2:2, 1:-1:2] = \
            (var_hkm_edge[1:-3:2, 1:-1:2] + var_hkm_edge[3:-1:2, 1:-1:2]) * 0.5
    var_hkm_edge[1:-1:1, 2:-2:2] = \
            (var_hkm_edge[1:-1:1, 1:-3:2] + var_hkm_edge[1:-1:1, 3:-1:2]) * 0.5

    # out edges
    var_hkm_edge[0,1:-1:1] = var_hkm_edge[1,1:-1:1] * 2.0 \
            - var_hkm_edge[2,1:-1:1]
    var_hkm_edge[-1,1:-1:1] = var_hkm_edge[-2,1:-1:1] * 2.0 \
            - var_hkm_edge[-3,1:-1:1]
    var_hkm_edge[:,0]  = var_hkm_edge[:,1]  * 2.0 - var_hkm_edge[:,2]
    var_hkm_edge[:,-1] = var_hkm_edge[:,-2] * 2.0 - var_hkm_edge[:,-3]

    # half km center
    var_hkm = (var_hkm_edge[:-1, :-1] + var_hkm_edge[1:, :-1] + 
               var_hkm_edge[:-1, 1:]  + var_hkm_edge[1:, 1:]) * 0.25

    return var_hkm, var_hkm_edge
#
#------------------------------------------------------------------------------
#
def modis_coastal_mask_pixel(all_ref):
    """ Implementation of Li et al. (2003)
    Li, Rong-Rong and Kaufman, Y. J. and Gao, Bo-Cai and Davis, C. O.,
    Remote sensing of suspended sediments and shallow coastal waters,
    Geoscience and Remote Sensing, IEEE Transactions on.

    Parameteres
    -----------
    all_ref : 1-D numpy array
        TOA reflectance at 0.466 um, 0.554 um, 1.241 nm,
        1.628 um, and 2.113 um.

    Returns
    -------
    out_dict : dict
        keys are ref_466, diff_554, mask

    """

    out_dict = {}
    out_dict['ref_466'] = all_ref[0]

    ref_466_thre = 0.25

    diff_554_thre = 0.01

    band = np.array( [0.466, 1.241, 1.628, 2.113] )
    out_dict['regress_band'] = band

    # if any vaule in all_ref is nan, mask is set as True
    if np.any(np.isnan(all_ref)):

        out_dict['diff_554'] = np.nan
        out_dict['mask'] = True

    else:

        # TOA reflectance at 0.466 um, 1.241 nm, 1.628 um, and 2.113 um.
        ref = np.hstack([np.array([all_ref[0]]), all_ref[2:]])

        # power law fit
        slope, intercept, r_value, p_value, std_err = \
                stats.linregress(np.log(band), np.log(ref))

        # ref_554 calculated from power law
        cal_ref_554 = np.log(0.544) * slope + intercept
        cal_ref_554 = np.exp(cal_ref_554)

        # cal_ref
        cal_ref = np.log(band) * slope + intercept
        cal_ref = np.exp(cal_ref)
        out_dict['regress_band_cal_ref'] = cal_ref

        # observational ref_554
        ref_554 = all_ref[1]

        # diff_554
        diff_554 = ref_554 - cal_ref_554
        out_dict['diff_554'] = diff_554

        # mask
        if (diff_554 > diff_554_thre) and (out_dict['ref_466'] < ref_466_thre):
            out_dict['mask'] = True
        else:
            out_dict['mask'] = False

    return out_dict
#
#------------------------------------------------------------------------------
#
