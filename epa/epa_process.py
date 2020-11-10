"""
Created on October 09, 2020

@author: Yi Wang
"""
import numpy as np

#
#------------------------------------------------------------------------------
#
def find_date_data(df, curr_date, n_obs_thre=12):
    """ find data from *df* for *curr_date*

    Parameters
    ----------
    df : pandas DataFrame
        Daily EPA data for one site.
    curr_date : str
        'YYYY-MM-DD'
    n_obs_thre : int
         At least *n_obs_thre* of hourly observations are
         used to calculate daily mean.

    returns
    -------
    out_dict : dict

    """

    out_dict = {}

    date  = np.array(df['Date Local'])
    n_obs = np.array(df['Observation Count'], dtype=int)
    val   = np.array(df['Arithmetic Mean'], dtype=float)

    # filter data
    flag = np.logical_and(date == curr_date, n_obs >= n_obs_thre)

    val = val[flag]
    if len(val) == 0:
        out_dict['value'] = np.nan
    else:
        out_dict['value'] = val[0]
        if out_dict['value'] < 0.0:
            out_dict['value'] = np.nan

    return out_dict 

#
#------------------------------------------------------------------------------
#
