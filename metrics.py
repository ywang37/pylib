"""
Created on September 25, 2019

@author: Yi Wang
"""

import numpy as np
import scipy.stats as ss

def NMSE_mean(obs, est):
    """ Normalized mean squared error.
    (normalized by mean(obs) * mean(est))

    Parameters
    ----------
    obs : numpy array
        Observations
    est : numpy array
        Estimations

    Returns
    -------
    result : float

    """

    result = np.mean((est-obs) ** 2) / (np.mean(est) * np.mean(obs))

    return result

def NMSE_var(obs, est):
    """ Normalized mean squared error.
    (normalized by variance(obs))

    Parameters
    ----------
    obs : numpy array
        Observations
    est : numpy array
        Estimations

    Returns
    -------
    result : float

    """

    result = np.mean((est-obs) ** 2) / np.var(obs, ddof=1)

    return result

def MSE(obs, est):
    """ Mean squared error.

    Parameters
    ----------
    obs : numpy array
        Observations
    est : numpy array
        Estimations

    Returns
    -------
    result : float

    """

    result = np.mean((est-obs) ** 2)

    return result

def RMSE(obs, est):
    """ Root mean squared error.

    Parameters
    ----------
    obs : numpy array
        Observations
    est : numpy array
        Estimations

    Returns
    -------
    result : float

    """

    result = np.sqrt(np.mean((est-obs) ** 2))

    return result

def MB(obs, est):
    """ Mean bias.

    Parameters
    ----------
    obs : numpy array
        Observations
    est : numpy array
        Estimations

    Returns
    -------
    result : float

    """

    result = np.mean(est) - np.mean(obs)

    return result
    
def NMB(obs, est):
    """ Normalized mean bias

    Parameters
    ----------
    obs : numpy array
        Observations
    est : numpy array
        Estimations

    Returns
    -------
    result : float

    """

    result = (np.mean(est) - np.mean(obs)) / np.mean(obs)

    return result

def NCRMSE(obs, est, norm=True):
    """ (Normalized) centralized root mean squared error

    Parameters
    ----------
    obs : numpy array
        Observations
    est : numpy array
        Estimations
    norm : logical
        Normalized by np.std(obs)

    Returns
    -------
    crmse or ncrmse : float
        

    """

    # linear correlation
    r, p_value = ss.pearsonr(obs, est)

    # standard deviation
    obs_std = np.std(obs, ddof=0)
    est_std = np.std(est, ddof=0)

    # centralized root mean squared error
    crmse = np.sqrt(obs_std ** 2 + est_std ** 2 \
            - 2 * obs_std * est_std * r)

    # Normalized centralized root mean squared error
    ncrmse = crmse / obs_std

    if norm:
        return ncrmse
    else:
        return crmse







    






