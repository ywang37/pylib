"""
Created on March 19, 2020

@author: Yi Wang
"""

import os

#
#------------------------------------------------------------------------------
#
def download_met(res, model, year, month, root_dir=None):
    """ Download meteorological data for GEOS-Chem
    (ywang, 03/19/2020)

    Parameters
    ----------
    res : str
        Resolution. 'GEOS_2x2.5', 'GEOS_0.25x0.3125_CH', for example.
    model : str
        'GEOS_FP', 'MERRA2'
    year : str
        '2018' for example.
    month : str
        '05' for example
    root_dir : str or None
        Parent directory to save data.
        If None, root_dir is '/Dedicated/jwang-data/GCDATA/'

    """

    # root URL
    url_root = 'http://geoschemdata.computecanada.ca/ExtData/'

    # URL
    url = url_root + res + '/' + model + '/' + year + '/' + \
            month + '/'
    #print(url)

    # root directory
    if root_dir is None:
        root_dir = '/Dedicated/jwang-data/GCDATA/'
    if root_dir[-1] != '/':
        root_dir = root_dir + '/'

    # data directory
    data_year_dir = root_dir + res + '/' + model + '/' + year + '/' 
    data_month_dir = data_year_dir + month + '/'

    # check year directory
    if not os.path.isdir(data_year_dir):
        os.system('mkdir -p ' + data_year_dir)
        os.system('chmod g+w ' + data_year_dir)

    # check month directory
    if not os.path.isdir(data_month_dir):
        os.system('mkdir -p ' + data_month_dir)
        os.system('chmod g+w ' + data_month_dir)
    else:
        print(data_month_dir + ' already exists, data is not downloaded.')
        return

    # download
    cmd = 'wget -r -np -nH -N -R "*.html*" --cut-dirs=5 ' + \
            '"' + url + '" -P "' + data_month_dir + '"'
    os.system(cmd)

    # access right
    os.system('chmod g+w ' + data_month_dir + '*')

#
#------------------------------------------------------------------------------
#
def link_for_soil_temp(res, model, year, month, root_dir=None):
    """ All meteorological data are lined except 'A1' data.
    New 'A1' data contain soil temperatuer.
    (ywang, 03/19/2020)

    Parameters
    ----------
    res : str
        Resolution. 'GEOS_2x2.5', 'GEOS_0.25x0.3125_CH', for example.
    model : str
        'GEOS_FP', 'MERRA2'
    year : str
        '2018' for example.
    month : str
        '05' for example
    root_dir : str or None
        Parent directory for meterorological data.
        If None, root_dir is '/Dedicated/jwang-data/GCDATA/'

    """

    if model == 'GEOS_FP':
        model_file = 'GEOSFP'

    # root directory
    if root_dir is None:
        root_dir = '/Dedicated/jwang-data/GCDATA/'
    if root_dir[-1] != '/':
        root_dir = root_dir + '/'

    # Original data directory
    ori_data_dir = root_dir + res + '/' + model + '/' + year + '/' \
            + month + '/'
    print('Original directory is ' + ori_data_dir)

    # New data directory
    new_data_dir = root_dir + res + '/' + model + '_soil_T/' + year + '/' \
            + month + '/'
    print('New directory is ' + new_data_dir)

    # check directory
    if not (os.path.isdir(ori_data_dir)  and os.path.isdir(new_data_dir)):
        print('At least one directory does not exist.')
        print('Link is not done.')

    # link
    os.system('ln -sf ' + ori_data_dir + model_file + '.' + year + month +
            '*.A3*.*.nc ' + new_data_dir)
    os.system('ln -sf ' + ori_data_dir + model_file + '.' + year + month +
            '*.I3.*.nc ' + new_data_dir)




#
#------------------------------------------------------------------------------
#
