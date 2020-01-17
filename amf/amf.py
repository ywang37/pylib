"""
Created on Janunary 16, 2020

@author: Yi Wang
"""

import numpy as np
from scipy import interpolate
from tqdm import tqdm

from mylib.model.trop import tropospheric_layer_column_one

#
#------------------------------------------------------------------------------
#
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
#
#------------------------------------------------------------------------------
#
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
#
#------------------------------------------------------------------------------
#
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

    # fit a function
    func = interpolate.interp1d(SW_AK_press, SW_AK)

    # Avoid extrapolation.
    tmp_new_press = np.array(new_press)
    tmp_new_press[tmp_new_press>SW_AK_press[0]]  = SW_AK_press[0]
    tmp_new_press[tmp_new_press<SW_AK_press[-1]] = SW_AK_press[-1]

    # Interpolation
    new_SW_AK = func(tmp_new_press)

    return new_SW_AK
#
#------------------------------------------------------------------------------
#
def AMF_trop_one(layer_val, PEdge_Bot, SW_AK, SW_AK_press, ind_l=None,
        P_tropopause=None, var='AMF'):
    """ Calculate tropospheric
    (1) Air mass factor if *var* is 'AMF', or
    (2) Vertical column density that can be compared to the 
        satellite retrieval if *var* is 'VCD_AK'

    Parameters
    ----------
    layer_val : 1-D array
        Species columnar cencentrations at every layer.
        Unit is like molec/cm^2, DU, ...
    PEdge_Bot : 1-D array or None
        Bottom edge pressure of every layer.
        Layer index starts from 0.
        Pressure decreases as layer index increases.
    SW_AK: 1-D array
        Scattering weight or averaging kernel
    SW_AK_press : 1-D array
        Pressure of *SW_AK*
    ind_l : float or None
        Layer index of tropospause.  Layer index starts from 0.
        If *ind_l* is None, *P_tropopause* should be given.
    P_tropopause : float or None
        Tropopause pressure. 
        If *P_tropopause* is None, *ind_l* should be given.
    var : str (default is 'AMF')
        'AMF': *SW_AK* is scattering weight, and air mass factor
            is calculated.
        'VCD_AK': *SW_AK* is averaging kernel, and vertical column 
            density that can be compared to the satellite retrieval.

    Returns
    -------
    out_dict : dict
        'AMF': Air mass factor if *var* is 'AMF'.
        'VCD_AK': Model tropospheric vertical column density that can 
            be compared to the satellite retrieval.
        'VCD' : Model tropospheric vertical column density 

    """

    out_dict = {}

    # get data of species in the troposphere 
    trop_data = tropospheric_layer_column_one(layer_val, ind_l=ind_l,
            PEdge_Bot=PEdge_Bot, P_tropopause=P_tropopause, pressure=True)
    out_dict['VCD'] = trop_data['trop_col']

    # Interpolate scattering weight or averaging kernel
    # to pressure of each tropospheric levels.
    trop_press = trop_data['trop_PCenter']
    trop_SW_AK = SW_AK_intep_one(SW_AK, SW_AK_press, trop_press)

    # Calculate air mass factor or vertical column density 
    # that can be compared to the satellite retrieval.
    trop_layer_val = trop_data['trop_layer']
    if (var == 'AMF'):

        # Air mass factor
        AMF_data = AMF_one(trop_layer_val, trop_SW_AK)
        AMF = AMF_data['AMF']
        out_dict['AMF'] = AMF

    elif (var == 'VCD_AK'):

        # Vertical column density that can be compared to
        # the satellite retrieval.
        VCD_AK = VCD_AK_one(trop_layer_val, trpo_SW_AK)
        out_dict['VCD_AK'] = VCD_AK

    else:

        print(' - AMF_trop_one: *var* is {}'.format(var))
        print('  *var*\'s vaule can be: ', ['AMF', 'VCD_AK'])
        print('  If *var* is \'AMF\', *SW_AK* is scatter weight. ' +
                'And air mass factor is calculated')
        print('  If *var* is \'VCD_AK\', *SW_AK* is averaging kernel. ' +
                'And vertical column density that can be compared to ' +
                'the satellite retrieval is calculated.')

    return out_dict
#
#------------------------------------------------------------------------------
#
def AMF_trop(layer_val_arr, PEdge_Bot_arr, SW_AK_arr, 
        SW_AK_press, ind_l_arr=None,
        P_tropopause_arr=None, var='AMF',
        flag=None):
    """

    flag : 
    """

    # Field dimension
    dim = layer_val_arr.shape
    dim = dim[0:-1]

    # total number of grid
    N_grid = 1
    for i in range(len(dim)):
        N_grid *= dim[i]

    # flag
    if flag is None:
        flag = np.full(dim, True)

    # output dict
    out_dict = {}
    out_dict['VCD'] = np.full(dim, np.nan)
    if (var in ['AMF', 'VCD_AK']):
        out_dict[var] = np.full(dim, np.nan)
    else:
        print(' - AMF_trop: *var* is {}'.format(var))
        print('  *var*\'s vaule can be: ', ['AMF', 'VCD_AK'])
        print('  If *var* is \'AMF\', *SW_AK_arr* is scatter weight. ' +
                'And air mass factor is calculated')
        print('  If *var* is \'VCD_AK\', *SW_AK_arr* is averaging kernel. ' +
                'And vertical column density that can be compared to ' +
                'the satellite retrieval is calculated.')


    # iterate every grid
    for i in tqdm(range(N_grid)):

        # get index
        if (len(dim) == 1):
            ind = i
        elif (len(dim) == 2):
            ind = (i // dim[1], i % dim[1])
        else:
            print(' - AMF_trop: dimension error.')
            exit()

        print(i, ind)

        # only process data with flag is True
        if (not flag[ind]):
            pass

        # prepare parameters for AMF_trop_one function.
        layer_val = layer_val_arr[ind,:]
        PEdge_Bot = PEdge_Bot_arr[ind,:]
        SW_AK     = SW_AK_arr[ind, :]
        if ind_l_arr is None:
            ind_l = None
        else: 
            ind_l = ind_l_arr[ind]
        if P_tropopause_arr is None:
            P_tropopause = None
        else:
            P_tropopause = P_tropopause_arr[ind]

        # calculate AMF or VCD_AK
        data_one = AMF_trop_one(layer_val, PEdge_Bot, SW_AK, SW_AK_press, 
                ind_l=ind_l, P_tropopause=P_tropopause, var=var)

        # save data
        out_dict[var][ind]   = data_one[var]
        out_dict['VCD'][ind] = data_one['VCD'] 

    return out_dict
#
#------------------------------------------------------------------------------
#
