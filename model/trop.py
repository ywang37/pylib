"""
Created on Janunary 15, 2020

@author: Yi Wang
"""

def tropospheric_column(layer_val, ind_l, flag=None):
    """ Calculate tropospheric column.

    Parameters
    ----------
    layer_val : 1-D, 2-D, or 3-D array
        Species columnar cencentrations at every layer.
        Unit is like molec/cm^2, DU, ...
        The last dimension is vertical layer.
    ind_l : array or scalar, dimension is one less than that of layer_val.
        Layer index of tropospause.  Layer index starts from 0.

    Returns
    -------
    trop_col : array or scalar, dimension is one less than that of layer_val.
        Tropospheric column

    """

    # index ceil
    ind_ceil = np.ceil(ind_l)
    ind_ceil = ind_ceil.astype(int)

    # ratio of tropospause layer in troposphere
    trop_pau_ratio = ind_l - (ind_ceil - 1.0)

    # one pixel
    if layer_val.ndim == 1:
        trop_col = np.sum(layer_val[0:ind_ceil]) + \
                (layer_val[ind_ceil] * trop_pau_ratio)

    # 1-D pixels
    if layer_val.ndim == 2:
        trop_col = np.zeros(layer_val[0:1])
        for i in range(trop_col.shape[0]):
            if flag[i]:
                ind_c = ind_ceil[i]
                trop_col[i] = np.sum(layer_val[i,0:ind_c]) + \
                        (layer_val[i,ind_c] * trop_pau_ratio[i])

    # 2-D pixels
    if layer_val.ndim == 2:
        trop_col = np.zeros(layer_val[0:2])
        for i in range(trop_col.shape[0]):
            for j in range(trop_col.shape[1]):
                if flag[i,j]:
                    ind_c = ind_ceil[i,j]
                    trop_col[i,j] = np.sum(layer_val[i,j,0:ind_c]) + \
                            (layer_val[i,j,ind_c] * trop_pau_ratio[i,j])

    return trop_col

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

    ind = i - 1.0 - ( (P_tropopause - PEdge_Bot[i]) / (PEdge_Bot[i-1] - PEdge_Bot[i]) )

    return ind
