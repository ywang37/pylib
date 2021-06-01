"""
Created on May 05, 2021

@author: Yi Wang
"""

import numpy as np
import pandas as pd

#
#------------------------------------------------------------------------------
#
def read_improve(filename, skiprows=16951):
    """ Read original IMPROVE data.
    (ywang, 05/13/2021)

    Parameters
    ----------
    filename : str
        IMPROVE filename
    skiprows : int or None
        Number of lines to skip (int) at the start of the file.

    Return
    df : DataFrame
        
    """

    df = pd.read_csv(filename, skiprows=skiprows)

    return df
#
#------------------------------------------------------------------------------
#
