"""
Created on Janunary 16, 2020

@author: Yi Wang
"""

import copy
import numpy as np
from scipy import interpolate
from tqdm import tqdm

from mylib.grid_utility import get_pressure_index
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
def correction_factor_one(AK, S_aprior, S_new):
    """ Calculate correction factor for satellite vertial column
    density when new shape facor is available.
    Here AK is not average kernel.
    I'm not sure if AK is scattering_weight / AMF_sat

    Parameters
    ----------
    AK : 1-D array
        Averageing kernel
    S_aprior : 1-D array
        aprior shape used in the trace gas retrieving procees
        if V is vertical column desnity, V * S_aprior is verical
        column density at every layer. np.sum(S_aprior) is 1.0
    S_new : 1-D array
        new shape facor

    Returns
    -------
    corrn_factor : float
        correction factor.
        If V_sat is satellite VCD retrieval, V_sat * corrn_factor is
        a new corrected satellite VCD retrieval.

    """

    # check AK
    if np.any( np.isnan(AK) ):
        print(' - correction_factor_one: nan exists in AK.')
        print(' AK is ', AK)
        exit()

    # check S_aprior
    if np.any( np.isnan(S_aprior) ):
        print(' - correction_factor_one: nan exists in S_aprior.')
        print(' S_aprior is ', S_aprior)
        exit()

    # check S_new
    if np.any( np.isnan(S_new) ):
        print(' - correction_factor_one: nan exists in S_new.')
        print(' S_new is ', S_new)
        exit()

    corrn_factor = 1.0 + np.sum( (AK - 1.0) * (S_aprior - S_new) )

    return corrn_factor
#
#------------------------------------------------------------------------------
#
def shape_factor_intep_one(layer_val_in, press_edge_in, new_press_edge):
    """ Calculate shape factor in *new_press_edge* 
    (Yi Wang, 01/29/2020)

    Parameters
    ----------
    layer_val_in : 1-D array
        Species columnar cencentrations at every layer.
        Unit is like molec/cm^2, DU, ...
    press_edge_in : 1-D array
        Pressure edge of *layer_val*
        Thus len(press_edge) == (len(layer_val) + 1) is usually True.
        Sometimes,  top level pressure is not provide, then 
        len(press_edge) == len(layer_val) is usually True.
        The first level is bottom.
    new_press_edge : 1-D array
        Shape factor is calculated at these pressure layers.
        The first level is bottom

    Returns
    -------
    out_dict : dict
        keys and values:
            'new_shape_factor' : 1-D array
                Shape factor at *new_press_edge*
            'new_layer_val' : 1-D array
                VCD at *new_press_edge*

    """

    out_dict = {}

    # copy data
    layer_val  = copy.deepcopy(layer_val_in)
    press_edge = copy.deepcopy(press_edge_in)

    # If the lowest level of *new_press_edge* has higher pressure
    # than the that of the lowest level of press_edge, we artificall
    # extend *layer_val* to the lowest level of *new_press_edge*
    press_bot     = press_edge[0]
    new_press_bot = new_press_edge[0]
    if (new_press_bot > press_bot):

        # extend pressure level
        press_edge = np.insert(press_edge, 0, new_press_bot)

        # linearly extend *layer_val* the lowest level
        press_diff   = press_edge[1] - press_edge[2]
        press_extend = press_edge[0] - press_edge[1]
        ratio = press_extend / press_diff
        val_extend = ratio * layer_val[0]
        layer_val = np.insert(layer_val, 0, val_extend)

    # Check if pressure levels are in descending order.
    # check press_edge
    for i in range(len(press_edge)-1):
        if (press_edge[i+1] >= press_edge[i]):
            print(' - shape_factor_intep: press_edge is not in ' + 
                    'descending order.')
            print('press_edge is: ')
            print(press_edge)
            exit()
    # check new_press_edge
    for i in range(len(new_press_edge)-1):
        if (new_press_edge[i+1] >= new_press_edge[i]):
            print(' - shape_factor_intep: new_press_edge is not in ' +
                    'descending order.')
            print('new_press_edge is: ')
            print(new_press_edge)
            exit()

    # Interpolate layer_val to new layers defined by new_press_edge,
    # and save the result to new_layer_val.
    new_layer_val = np.full((len(new_press_edge)-1,), np.nan)
    for i in range(len(new_layer_val)):

        # pressure bottom and top level for each layer
        press_floor   = new_press_edge[i]
        press_ceiling = new_press_edge[i+1]

        # layer index in layer_val 
        j_floor   = get_pressure_index(press_edge, press_floor  )
        j_ceiling = get_pressure_index(press_edge, press_ceiling)

        # calcualte VCD within [press_floor, press_ceiling]
        if (j_floor == j_ceiling):
        # New layer is subset of an old layer

            new_layer_val[i] = layer_val[j_floor] * \
                    ( (press_floor         - press_ceiling        ) / \
                      (press_edge[j_floor] - press_edge[j_floor+1])     )

        elif (j_floor < j_ceiling):
        # New layer is subset of at least two old layers

            # bottom part of the new layer
            bot_part = layer_val[j_floor] * \
                    ( (press_floor         - press_edge[j_floor+1]) / \
                      (press_edge[j_floor] - press_edge[j_floor+1])     )

            # top part of the new layer
            top_part = layer_val[j_ceiling] * \
                    (  (press_edge[j_ceiling] - press_ceiling          ) / \
                       (press_edge[j_ceiling] - press_edge[j_ceiling+1])    )

            # add bottom part and top part, and assign the sum
            # to new layer
            new_layer_val[i] = bot_part + top_part

            # add center part of the new layer, if exists
            if ( (j_ceiling - j_floor) > 1 ):
                new_layer_val[i] += np.sum( layer_val[j_floor+1:j_ceiling] )

        else:
        # error
            print(' - shape_factor_intep: j_floor is {}, which is larger' + \
                    ' than j_ceiling {}'.format(j_floor, j_ceiling))

    # VCD at new_press_edge
    out_dict['new_layer_val'] = new_layer_val

    # new shape factor at new_press_edge
    out_dict['new_shape_factor'] = new_layer_val / np.sum(new_layer_val)

    return out_dict
#
#------------------------------------------------------------------------------
#
def shape_factor_correction_factor_one(layer_val, press_edge, 
        new_press_edge_in, AK_in, S_aprior_in, nan_flag=True):
    """ Call shape_factor_intep_one function and correction_factor_one
    function to caulculate correction_factor for one pixel.
    (Yi Wang, 01/30/2020)

    Parameters
    ----------
    layer_val : 1-D array
        Species columnar cencentrations at every layer.
        Unit is like molec/cm^2, DU, ...
    press_edge : 1-D array
        Pressure edge of *layer_val*
        Thus len(press_edge) == (len(layer_val) + 1) is usually True.
        Sometimes,  top level pressure is not provide, then 
        len(press_edge) == len(layer_val) is usually True.
        The first level is bottom.
    new_press_edge_in : 1-D array
        Shape factor is calculated at these pressure layers.
        The first level is bottom
    AK_in : 1-D array
        Averageing kernel in the *new_press_edge* layers
    S_aprior_in : 1-D array
        aprior shape used in the trace gas retrieving procees
        if V is vertical column desnity, V * S_aprior is verical
        column density at every layer. np.sum(S_aprior) is 1.0
        *S_aprior* is in the *new_press_edge* layers
    nan_flag : bool (default is Ture)
        If nan_flag is True, np.nan in the low layers of new_press_edge,
        AK, and S_aprior will be removed.

    Returns
    -------
    out_dict : dict
        keys and values:
            'new_shape_factor' : 1-D array
                Shape factor at *new_press_edge*
            'new_layer_val' : 1-D array
                VCD at *new_press_edge*
            'corrn_factor' : float
                Correction factor

    """

    # remove nan
    if nan_flag:
        # find index of first valid element
        for i in range(len(new_press_edge_in)):
            if not np.isnan(new_press_edge_in[i]):
                break
        # new_press_edge
        new_press_edge = copy.copy(new_press_edge_in[i:])
        if np.any( np.isnan(new_press_edge) ):
            print(' - shape_factor_correction_factor_one: ' + 
                    'nan exists in new_press_edge.')
            print(' new_press_edge_in is ', new_press_edge_in)
            print(' new_press_edge is ', new_press_edge)
            exit()
        # AK
        AK = copy.copy(AK_in[i:])
        if np.any( np.isnan(AK) ):
            print(' - shape_factor_correction_factor_one: ' + 
                    'nan exists in AK.')
            print(' AK_in is ', AK_in)
            print(' AK is ', AK)
            exit()
        # S_aprior
        S_aprior = copy.copy(S_aprior_in[i:])
        if np.any( np.isnan(S_aprior) ):
            print(' - shape_factor_correction_factor_one: ' + 
                    'nan exists in S_aprior.')
            print(' S_aprior_in is ', S_aprior_in)
            print(' S_aprior is ', S_aprior)
            exit()
    else:
        new_press_edge = copy.copy(new_press_edge_in)
        AK             = copy.copy(AK_in)
        S_aprior       = copy.copy(S_aprior_in)

    # calculate new shape factor in the *new_press_edge* layer
    shape_dict = shape_factor_intep_one(layer_val, press_edge, new_press_edge)

    # calculate correction factor
    new_shape_factor = shape_dict['new_shape_factor']
    corrn_factor = correction_factor_one(AK, S_aprior, new_shape_factor)

    # output data
    out_dict = {}
    out_dict['new_shape_factor'] = new_shape_factor
    out_dict['new_layer_val']    = shape_dict['new_layer_val']
    out_dict['corrn_factor']     = corrn_factor

#    #if np.sum(layer_val)>0.2:
#    #if np.sum(layer_val[0])>0.06:
#    if (np.sum(np.isnan(new_press_edge_in))>3):
#        print('layer_val')
#        print(layer_val)
#        print('press_edge')
#        print(press_edge)
#        print('new_press_edge_in')
#        print(new_press_edge_in)
#        print('AK_in')
#        print(AK_in)
#        print('S_aprior_in')
#        print(S_aprior_in)
#        print('new_press_edge')
#        print(new_press_edge)
#        print('AK')
#        print(AK)
#        print('S_aprior')
#        print(S_aprior)
#        print('new_shape_factor')
#        print(new_shape_factor)
#        print('new_layer_val')
#        print(out_dict['new_layer_val'])
#        print('corrn_factor={}'.format(corrn_factor))
#        exit()

    return out_dict
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
    for one pixel

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
#    print('AMF_trop_one')
#    print(P_tropopause)
#    print(PEdge_Bot)
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
    """ Calculate tropospheric
    (1) Air mass factor if *var* is 'AMF', or
    (2) Vertical column density that can be compared to the 
        satellite retrieval if *var* is 'VCD_AK'
    for 1-D or 2-D field.

    Parameters
    ----------
    layer_val_arr : 2-D or 3-D array
        The elements are layer_val parameter for AMF_trop_one function.
    PEdge_Bot_arr : 2-D or 3-D array or None
        The elements are PEdge_Bot parameter for AMF_trop_one function.
    SW_AK_arr : 2-D or 3-D array
        The elements are SW_AK parameter for AMF_trop_one function.
    SW_AK_press : 1-D array
        Pressure of *SW_AK_arr*
    ind_l_arr : 1-D or 2-D array or None
        The elements are ind_l parameter for AMF_trop_one function.
    P_tropopause_arr : 1-D or 2-D array or None
        The elements are P_tropopause parameter for AMF_trop_one function.
    var : str (default is 'AMF')
        'AMF': *SW_AK* is scattering weight, and air mass factor
            is calculated.
        'VCD_AK': *SW_AK* is averaging kernel, and vertical column 
            density that can be compared to the satellite retrieval.
    flag : 1-D or 2-D bool array or None
        Only precess then pixel with True flag. If flag is None,
        all pixels are processed.

    Returns
    -------
    out_dict : dict
        'AMF': Air mass factor if *var* is 'AMF'.
        'VCD_AK': Model tropospheric vertical column density that can 
            be compared to the satellite retrieval.
        'VCD' : Model tropospheric vertical column density 

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
    #for i in tqdm(range(N_grid)):
    for i in range(N_grid):

        # get index
        if (len(dim) == 1):
            ind = i
        elif (len(dim) == 2):
            ind = (i // dim[1], i % dim[1])
        else:
            print(' - AMF_trop: dimension error. 1')
            print('dim is: ')
            print(dim)
            exit()

        # only process data with flag is True
        if  (not flag[ind]):
            continue


        # prepare parameters for AMF_trop_one function.
        if (len(dim) == 1):
            layer_val = layer_val_arr[ind,:]
            PEdge_Bot = PEdge_Bot_arr[ind,:]
            SW_AK     = SW_AK_arr[ind, :]
        elif (len(dim) == 2):
            layer_val = layer_val_arr[ind[0],ind[1],:]
            PEdge_Bot = PEdge_Bot_arr[ind[0],ind[1],:]
            SW_AK     = SW_AK_arr[ind[0],ind[1], :]
        else:
            print(' - AMF_trop: dimension error. 2')
            exit()
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
def shape_factor_correction_factor(layer_val_arr, press_edge_arr,
        new_press_edge_arr, AK_arr, S_aprior_arr, nan_flag=True,
        flag=None):
    """ Call shape_factor_correction_factor_one function to calculate
    correction factors for 1-D or 2-D field.
    (Yi Wang, 01/30/2020)

    Parameters
    ----------
    layer_val_arr : 2-D or 3-D array
        The elements are layer_val parameter for 
        shape_factor_correction_factor_one function.
    press_edge_arr : 2-D or 3-D array
        The elements are press_edge parameter for
        shape_factor_correction_factor_one function.
    new_press_edge_arr : 2-D or 3-D array
        The elements are new_press_edge_in parameter for
        shape_factor_correction_factor_one function.
    AK_arr : 2-D or 3-D array
        The elements are AK_in parameter for
        shape_factor_correction_factor_one function.
    S_aprior_arr : 2-D or 3-D array
        The elements are S_aprior_in parameter for
        shape_factor_correction_factor_one function.
    nan_flag : bool (default is Ture)
        passed to shape_factor_correction_factor_one function.
        If nan_flag is True, np.nan in the low layers of 
        new_press_edge_arr, AK_arr, and S_aprior_arr will be removed.
    flag : 1-D or 2-D bool array or None
        Only precess then pixel with True flag. If flag is None,
        all pixels are processed.

    Returns
    -------
    out_dict : dict

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
    out_dict['corrn_factor'] = np.full(dim, np.nan)

    # iterate every grid
    #for i in tqdm(range(N_grid)):
    for i in range(N_grid):

        # get index
        if (len(dim) == 1):
            ind = i
        elif (len(dim) == 2):
            ind = (i // dim[1], i % dim[1])
        else:
            print(' - shape_factor_correction_factor: dimension error. 1')
            exit()

        # only process data with flag is True
        if  (not flag[ind]):
            continue

        # prepare parameters for AMF_trop_one function.
        if (len(dim) == 1):
            layer_val      = layer_val_arr[ind,:]
            press_edge     = press_edge_arr[ind,:]
            new_press_edge = new_press_edge_arr[ind,:]
            AK             = AK_arr[ind,:]
            S_aprior       = S_aprior_arr[ind,:]
        elif (len(dim) == 2):
            layer_val      = layer_val_arr[ind[0],ind[1],:]
            press_edge     = press_edge_arr[ind[0],ind[1],:]
            new_press_edge = new_press_edge_arr[ind[0],ind[1],:]
            AK             = AK_arr[ind[0],ind[1],:]
            S_aprior       = S_aprior_arr[ind[0],ind[1],:]
        else:
            print(' - shape_factor_correction_factor: dimension error. 2')
            exit()

        # calculate correction factor
        data_one = shape_factor_correction_factor_one(layer_val, press_edge,
                new_press_edge, AK, S_aprior, nan_flag=nan_flag)

        # save data
        out_dict['corrn_factor'][ind] = data_one['corrn_factor'] 

    return out_dict
#
#------------------------------------------------------------------------------
#
def test_driver_shape_factor_intep_one(old_layer_val, old_press_edge, 
        new_press_edge):

    out_dict = shape_factor_intep_one(old_layer_val, old_press_edge,
            new_press_edge)

    print('old_press_edge: ')
    print(old_press_edge)
    print('old_layer_val: ')
    print(old_layer_val)
    print('sum of old_layer_val = {}'.format(np.sum(old_layer_val)))
    print('new_press_edge: ')
    print(new_press_edge)
    print('new_layer_val: ')
    print(out_dict['new_layer_val'])
    print('sum of new_layer_val = ' + \
            '{}'.format(np.sum(out_dict['new_layer_val'])))
    print('new_shape_factor: ')
    print(out_dict['new_shape_factor'])
    print('sum of new_shape_factor = ' + \
            '{}'.format(np.sum(out_dict['new_shape_factor'])))
#
#------------------------------------------------------------------------------
#
if ('__main__' == __name__):

    n_dash = 30

    ###############################
    # Test shape_factor_intep_one
    ###############################

    print('-'*n_dash + 'Start testing shape_factor_intep_one' + '-'*n_dash)

    # old pressure edges
    old_press_edge = np.array(
            [1013.25, 1010.0, 1000.0, 975.0, 945.0, 900.0, 825.0, 700.0,
                550.0, 200.0, 100.0]
            )

    # old layer values
    old_layer_val = np.array(
            [1.2, 1.8, 2.6, 2.2, 1.9, 1.0, 0.6, 0.5, 0.2, 0.1]
            )

    # new pressure edges 1
    print('--- test 1 ---')
    new_press_edge_1 = np.array(old_press_edge)
    test_driver_shape_factor_intep_one(old_layer_val, old_press_edge,
            new_press_edge_1)

    # new pressure edges 2
    print('--- test 2 ---')
    new_press_edge_2 = np.array(
            [1013.25, 1010.0, 980.0, 800.0, 100.0]
            )
    test_driver_shape_factor_intep_one(old_layer_val, old_press_edge,
            new_press_edge_2)

    # new pressure edges 3
    print('--- test 3 ---')
    new_press_edge_3 = np.array(
            [1013.0, 1011.0, 920.0, 800.0, 700.0, 200.0, 100.0]
            )
    test_driver_shape_factor_intep_one(old_layer_val, old_press_edge,
            new_press_edge_3)

    # new pressure edges 4
    print('--- test 4 ---')
    new_press_edge_4 = np.array(
            [1015.0, 1011.0, 920.0, 800.0, 700.0, 560.0, 400.0, 300.0]
            )   
    test_driver_shape_factor_intep_one(old_layer_val, old_press_edge,
            new_press_edge_4)

    # new pressure edges 5
    print('--- test 5 ---')
    new_press_edge_5 = np.array(
            [1013.0, 1011.0, 920.0, 800.0, 700.0, 560.0, 400.0, 300.0]
            )
    test_driver_shape_factor_intep_one(old_layer_val, old_press_edge,
            new_press_edge_5)

    print('-'*n_dash + 'End testing shape_factor_intep_one' + '-'*n_dash)
