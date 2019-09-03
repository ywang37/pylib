"""
Created on September 2, 2019

@author: Yi Wang
"""

import numpy as np

def cal_ae(band1, band2, aod1, aod2):
    """ Calculate angstrom exponent
    """

    ae = np.log(aod1 / aod2) / np.log(band2 / band1)

    return ae

def interpolate_aod(band1, aod1, ae, band2):
    """ Interpolate *aod1* at *band1* to *band2*
    """
    aod2 = aod1 * ((band1/band2) ** ae)

    return aod2
