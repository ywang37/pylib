"""
Created on September 7, 2019

@author: Yi Wang
"""

import datetime

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
