"""
Created on Janunary 7, 2020

@author: Yi Wang
"""

####################
# chemical constants
####################

# Avogadro constatn
avo = 6.02214076e23

# molecular weight
mol_N = 14.0 # Nitrigon
mol_O = 16.0 # Oxygen

mol_dryair = 28.97 # dry air

mol_NO = 30.0 # Nitric oxide

##################
# unit conversions
##################

# mol/cm^2 => molec/cm^2
mol_m2_to_molec_cm2 = avo / 1e4

# molec/cm^2 => DU
molec_cm2_to_DU = 1.0 / 2.69e16

# molec (NOx) => kg (N)
molec_to_kgN = 1.0 / avo * mol_N / 1e3

# kg NO => ng N
kg_NO_to_ng_N = 1e12 / mol_NO * mol_N
