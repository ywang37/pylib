"""
Created on Janunary 15, 2020

@author: Yi Wang
"""

from mylib.constants import molec_cm2_to_DU

#################
# unit conversion
#################

def vmr_to_molec_cm2(val_old, air_den, box_hgt, vmr=None):
    """ Volume mixing ratio to molec/cm^2

    Parameters
    ----------
    val_old :
        unit is ppbv, ppmv, or v/v, depends on *vmr*
    air_den :
        air density. unit is molce/cm^2
    box_hgt :
        box height. unit is m.
    vmr :
        vaule is 'ppbv', 'ppmv', 'v/v', or None.
        None means 'ppbv'. *val_old* unit should
        be consistent wiht *vmr*

    Returns
    -------
    val_new:
        unit is molec/cm^2

    """


    if vmr is None:
        vmr = 'ppbv'

    #
    if vmr == 'ppbv':
        scale = 1e-9
    elif vmr == 'ppmv':
        scale = 1e-6
    elif vmr == 'v/v':
        scale = 1.0
    else:
        print(' - vmr_to_molec_cm2: vmr value is ', vmr,
                '. it only can be one of ', 
                ['ppbv', 'ppmv', 'v/v'])

    # conversion
    val_new = val_old * scale * air_den * (box_hgt * 100.0)

    return val_new

def vmr_to_DU(val_old, air_den, box_hgt, vmr=None):
    """ Volume mixing ratio to DU

    Parameters
    ----------
    val_old :
        unit is ppbv, ppmv, or v/v, depends on *vmr*
    air_den :
        air density. unit is molce/cm^2
    box_hgt :
        box height. unit is m.
    vmr :
        vaule is 'ppbv', 'ppmv', 'v/v', or None.
        None means 'ppbv'. *val_old* unit should
        be consistent wiht *vmr*

    Returns
    -------
    val_new:
        unit is DU

    """

    # conversion
    val_new = vmr_to_molec_cm2(val_old, air_den, box_hgt, vmr=None) \
            * molec_cm2_to_DU

    return val_new
