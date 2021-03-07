"""
Created on September 7, 2019

@author: Yi Wang
"""
from astropy.time import Time
import datetime

#
#------------------------------------------------------------------------------
#
def get_day_of_year(curr_date):
    """ Get day of year of *curr_date*

    Parameters
    ----------
    curr_date : str
        'YYYY-MM-DD'

    Returns
    -------
    day_of_year : int
    
    """

    start_date = curr_date[0:4] + '-01-01'

    curr_date_d  = datetime.datetime.strptime(curr_date,  '%Y-%m-%d')
    start_date_d = datetime.datetime.strptime(start_date, '%Y-%m-%d')

    diff = curr_date_d - start_date_d

    day_of_year = diff.days + 1

    return day_of_year
#
#------------------------------------------------------------------------------
#
def day_of_year_to_date(year_doy):
    """ Given day of year, output date

    Parameters
    ----------
    year_doy : str
        'YYYYDOY', for example, '2015006', '2015126'

    returns
    -------
    curr_date : str
        'YYYY-MM-DD'
    """

    year = year_doy[0:4]
    doy  = int(year_doy[4:])

    start_date_d = \
            datetime.datetime.strptime(year + '-01-01', '%Y-%m-%d')

    curr_date_d = start_date_d + datetime.timedelta(days=(doy-1))

    curr_date = str(curr_date_d)[0:10]

    return curr_date
#
#------------------------------------------------------------------------------
#
def tai_st(ymdhms, st='1993-01-01T00:00:00'):
    """ TAI seconds since *st*

    Parameters
    ----------
    ymdhms : str
        'yyyy-mm-ddThh:MM:ss', for example '2021-03-04T11:12:00'
    st : str
        'yyyy-mm-ddThh:MM:ss', for example '1993-01-01T00:00:00'

    Returns
    -------
    TAI : float
        TAI since *st*

    """

    one_day_TAI = 24 * 3600E0
    TAI = (Time(ymdhms).jd - Time(st).jd) * one_day_TAI

    return TAI
#
#------------------------------------------------------------------------------
#
