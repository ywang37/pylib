"""
Created on Janunary 16, 2020

@author: Yi Wang
"""

import numpy as np
from scipy import interpolate

def AMF_one(layer_val, SW):
    """ Calculate air mass factor

    Parameters
    ----------
    layer_val : 1-D array
        Species columnar cencentrations at every layer.
        Unit is like molec/cm^2, DU, ...
    SW : 1-D array
        Scattering weight

    Returns
    -------
    out_dict : dict
        'SCD': Slant column density
        'AMF': Air mass factor
    """

    out_dict = {}

    # slant column density
    SCD = np.sum( layer_val * SW )
    out_dict['SCD'] = SCD

    # vertical column density
    VCD = np.sum( layer_val )

    # air mass factor
    AMF = SCD / VCD
    out_dict['AMF'] = AMF

    return out_dict

def VCD_AK_one(layer_val, AK):
    """ Calculate a vertical column density that can be compared
    to the satellite retrieval.

    Parameters
    ----------
    layer_val : 1-D array
        Species columnar cencentrations at every layer.
        Unit is like molec/cm^2, DU, ...
    AK : 1-D array
        Averaging kernel

    Returns
    -------
    V_AK : float
        a vertical column density that can be compared
        to the satellite retrieval.

    """

    V_AK = np.sum( layer_val * AK )

    return V_AK

def SW_AK_intep_one(SW_AK, SW_AK_press, new_press):
    """ Interpolate scattering weight or averaging kernel
    to new pressure levels.

    Parameters
    ----------
    SW_AK: 1-D array
        Scattering weight or averaging kernel
    SW_AK_press : 1-D array
        Pressure of *SW_AK*
    new_press : 1-D array
        New pressure levels where *SW_AK* is to be
        interpolated.

    Returns
    -------
    new_SW_AK : 1-D array
        Scattering weight or averaging kernel at *new_press* 
        pressure levels.

    """

    func = interpolate.interp1d(SW_AK_press, SW_AK)

    new_SW_AK = func(new_press)

    return new_SW_AK
