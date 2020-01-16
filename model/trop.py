"""
Created on Janunary 15, 2020

@author: Yi Wang
"""
#
#------------------------------------------------------------------------------
#
def tropospheric_layer_column_one(layer_val, ind_l=None, 
        PEdge_Bot=None, P_tropopause=None, pressure=True):
    """ Calculate columnar concentations of each layer in
    the troposphere and total tropospheric column. If both
    *PEdge_Bot* and *P_tropopause* are given, center and 
    edge pressures of each layer in the the troposphere are
    calculated.

    Parameters
    ----------
    layer_val : 1-D array
        Species columnar cencentrations at every layer.
        Unit is like molec/cm^2, DU, ...
    ind_l : float or None
        Layer index of tropospause.  Layer index starts from 0.
        If *ind_l* is None, *PEdge_Bot* and *P_tropopause* should
        be given.
    PEdge_Bot : 1-D array or None
        Bottom edge pressure of every layer.
        Layer index starts from 0.
        Pressure decreases as layer index increases.
    P_tropopause : float or None
        Tropopause pressure. 
        If *P_tropopause* is None, *PEdge_Bot* and *ind_l* should
        be given
    pressure : bool (default is True)
        Whether or not return edge and center pressures of each
        layer in the troposphere

    Returns
    -------
    out_dict : dict
        'trop_layer'   : tropospherc column
        'trop_col'     : columnar concentations of each layer
        If *pressure* is True.
        'trop_PEdge'   : Edge pressures of each layer in the troposphere.
        'trop_PCenter' : Center pressures of each layer in the troposphere.

    """

    out_dict = {}

    # If *ind_l* is None. calcualte it.
    # *PEdge_Bot* and *P_tropopause* should be given.
    if ind_l is None:
        if ( (PEdge_Bot is not None) and (P_tropopause is not None) ):
            ind_l = find_tropopause_index_one(PEdge_Bot, P_tropopause)
        else:
            print(' - tropospheric_layer_column_one: *ind_l* is {}. ' + 
                    'both *PEdge_Bot* and *P_tropopause* ' + 
                    'should be given.'.format(ind_l))
            exit()

    # index ceil
    ind_ceil = np.ceil(ind_l)
    ind_ceil = ind_ceil.astype(int)

    # ratio of tropospause layer in troposphere
    trop_pau_ratio = ind_l - (ind_ceil - 1.0)

    # columnar concentations of each layer
    trop_layer = np.append( layer_val[0:ind_ceil], 
            layer_val[ind_ceil] * trop_pau_ratio )
    out_dict['trop_layer'] = trop_layer

    # tropospherc column
    trop_col = np.sum(trop_layer)
    out_dict['trop_col'] = trop_col

    # center and edge pressures of each layer in the the troposphere.
    if pressure:

        # Tropopause
        if P_tropopause is None:

            if ( (PEdge_Bot is not None) and (ind_l is not None) ):
                #P_tropopause = \
                #        find_tropopause_pressure_one(PEdge_Bot, ind_l)
                print(' - tropospheric_layer_column_one: ' +
                'function find_tropopause_pressure_one ' + 
                'is not implemented yet.')
                exit()
            else:
                print(' - tropospheric_layer_column_one: ')
                print('   *pressure* is {}'.format(pressure))
                print('   *P_tropopause* is {}'.format(P_tropopause))
                print('   Both *PEdge_Bot* and *ind_l* should be given.')
                exit()
    
        trop_PEdge   = np.append( PEdge_Bot[0:ind_ceil+1], P_tropopause )
        trop_PCenter = (trop_PEdge[0:-1] + trop_PEdge[1:]) * 0.5
        out_dict['trop_PEdge']   = trop_PEdge
        out_dict['trop_PCenter'] = trop_PCenter

    return out_dict
#
#------------------------------------------------------------------------------
#
def find_tropopause_index_one(PEdge_Bot, P_tropopause):
    """ Find tropopause index for one pixel

    Parameters
    ----------
    PEdge_Bot : 1-D array
        Bottom edge pressure of every layer.
        Layer index starts from 0.
        Pressure decreases as layer index increases.
    P_tropopause : float
        Tropopause pressure.

    Returns
    -------
    ind : float
        Tropopause index.

    """

    if (P_tropopause > PEdge_Bot[0]) or (P_tropopause < PEdge_Bot[-1]):
        print(' - find_tropopause_index_one: !!! ERROR !!! ' + \
                'Tropopause pressure is out of range')
        print(' Tropopause pressure is {}'.format(P_tropopause))
        print(' Bottom edge pressure of every layer: ')
        print(PEdge_Bot)

    for i in range(len(PEdge_Bot)):
        if (P_tropopause >= PEdge_Bot[i]):
            break

    ind = i - 1.0 \
            - ( (P_tropopause - PEdge_Bot[i]) \
            / (PEdge_Bot[i-1] - PEdge_Bot[i]) )

    return ind
#
#------------------------------------------------------------------------------
#
