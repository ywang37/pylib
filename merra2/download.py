"""
Created on May 05, 2020

@author: Yi Wang
"""

import datetime
import os


#
#------------------------------------------------------------------------------
#
def tavg1_url(yyyymmdd, collections=None):
    """ Generate a list that contains the URLs of MERRA-2
    one-day asm tavg1 data if collections is None.
    But users can specific variables now.
    (ywang, 05/05/2020)

    Parameters
    ----------
    yyyymmdd : str
        Date, '20180501' for example
    collections : list
        list of variables.  

    Returns
    -------
    urls : list
        All URLs
    """

    # URL of GEOS-FP data
    url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/'

    #
    if collections is None:
        collections = ['tavg1_2d_flx_Nx', 'tavg1_2d_lnd_Nx',
                       'tavg1_2d_rad_Nx', 'tavg1_2d_slv_Nx']

    # collecation directory
    coln_dir = {}
    coln_dir['tavg1_2d_flx_Nx'] = 'M2T1NXFLX.5.12.4'
    coln_dir['tavg1_2d_lnd_Nx'] = 'M2T1NXLND.5.12.4'
    coln_dir['tavg1_2d_rad_Nx'] = 'M2T1NXRAD.5.12.4'
    coln_dir['tavg1_2d_slv_Nx'] = 'M2T1NXSLV.5.12.4'

    #################
    # Generate URLs
    #################

    urls = []

    # date
    yyyy = yyyymmdd[0:4]
    mm   = yyyymmdd[4:6]
    dd   = yyyymmdd[6:8]

    # stream number
    yyyy_int = int(yyyy)
    if (yyyy_int >= 1980) and (yyyy_int <=1991):
        stream_num = '100'
    elif (yyyy_int >= 1992) and (yyyy_int <=2000):
        stream_num = '200'
    elif (yyyy_int >= 2001) and (yyyy_int <=2010):
        stream_num = '300'
    elif (yyyy_int >= 2011):
        stream_num = '400'

    for coln in collections:

        filename = url + coln_dir[coln] + '/' + yyyy + '/' + mm + \
                '/MERRA2_' + stream_num  + '.' + coln  + '.' + \
                yyyymmdd  + '.nc4'
        urls.append(filename)

    return urls
#
#------------------------------------------------------------------------------
#
def get_tavg1(yyyymmdd, data_dir, retry=10,
        collections=None, 
        user='yi-wang-4@uiowa.edu', 
        password=''):
    """ Download MERRA-2 one-day asm tavg1 data, if collections is None.
    But users can specific variables now.
    (ywang, 05/05/2020)

    Parameters
    ----------
    yyyymmdd : str
        Date, '20180501' for example
    data_dir : str
        The directory where all data will be saved to.
    retry : int
        Maximum number of retry if a file is not downloaded
        successfully. 0 means it will not check whether or not
        success.
    collections : list
        list of variables. 
    user : str
        MERRA-2 website user name.
    password : str
        MERRA-2 website password.

    Returns
    -------
    out_dict : dict

    """

    out_dict = {}

    # Save filename to be downloaded to tmp_tavg1_files.txt
    tmp_file = './tmp_tavg1_files.txt'
    f = open(tmp_file, 'w')
    urls = tavg1_url(yyyymmdd, collections=collections)
    out_dict['urls'] = urls
    for filename in urls:
        f.write(filename + '\n')

    # Close tmp_tavg1_files.txt
    f.close()

    # Download data
    cmd = 'wget -i ' + tmp_file + ' -P ' + data_dir + \
            ' --user=' + user + \
            ' --password=' + password
    os.system(cmd)

    # Remove tmp_tavg1_files.txt
    os.system('rm -rf ' + tmp_file)

    # check
    check_tavg1(urls, data_dir, retry=retry, user=user, password=password)

    return out_dict
#
#------------------------------------------------------------------------------
#
def check_tavg1(yyyymmdd, data_dir, retry=10, collections=None,
        user='yi-wang-4@uiowa.edu',
        password=''):
    """ Check if any MERRA-2 one-day asm tavg1 data (collections is None)
    are not downloaded successfully. If so, download them again.
    But users can specific variables now.
    (ywang, 05/05/2020)

    Parameters
    ----------
    yyyymmdd : str or list
        If str, Date, '20180501' for example
        If list, it contains the URLs of GEOS-FP 
        one-day asm tavg1 data.
    data_dir : str
        The directory where all data will be saved to.
    collections : list
            list of variables.
    user : str
        MERRA-2 website user name.
    password : str
        MERRA-2 website password.

    Returns
    -------

    """

    # make sure ncdump command is available
    tmp0 = os.system('ncdump > /dev/null 2>&1')
    if tmp0 != 0:
        print(' - check_tavg1: ncdump command is not available.')
        exit()

    if isinstance(yyyymmdd, list):
        urls = yyyymmdd
    elif isinstance(yyyymmdd, str):
        urls = tavg1_url(yyyymmdd, collections=collections)
    else:
        print(' - check_tavg1: yyyymmdd parameter error.')
        exit()

    if data_dir[-1] != '/':
        data_dir = data_dir + '/'

    # check every file.
    for url in urls:

        filename = data_dir + url.split('/')[-1] 
        nc_cmd = 'ncdump -h ' + filename + ' > /dev/null 2>&1'

        itry = 0
        while itry < retry:

            itry += 1

            # check a file
            tmp1 = os.system(nc_cmd)
            if tmp1 != 0:
                wget_cmd = 'wget ' + url + ' -O ' + filename + \
                        ' --user=' + user + \
                        ' --password=' + password
                print('Retry "' + wget_cmd + '" {:} time(s).'.format(itry))
                os.system(wget_cmd)

        # Fail after maxium retry.
        if (os.system(nc_cmd) != 0) and (itry == retry):
            print(url + ' failed after {} tries'.format(itry))
#
#------------------------------------------------------------------------------
#
def get_tavg1_month(yyyymm, root_data_dir,
        collections=None,
        user='yi-wang-4@uiowa.edu',
        password=''):
    """ Download MERRA2 one-month data.
    (ywang, 05/05/2020)

    Parameters
    ----------
    yyyymm : str
        Month, '201805' for example.
    root_data_dir : str
        The root directory where all data will be saved to.
    collections : list
        list of variables.
    user : str
        MERRA-2 website user name.
    password : str
        MERRA-2 website password.

    """

    if root_data_dir[-1] != '/':
        root_data_dir = root_data_dir + '/'

    currDate_D = datetime.datetime.strptime(yyyymm + '01', '%Y%m%d')

    while True:

        # Date
        currDate = str(currDate_D)
        yyyy = currDate[0:4]
        mm   = currDate[5:7]
        dd   = currDate[8:10]
        yyyymmdd = yyyy + mm + dd

        # Directory
        data_dir = root_data_dir + yyyy + '/' + mm
        if not os.path.isdir(data_dir):
            os.system('mkdir -p ' + data_dir)
        os.system('rm -f ' + data_dir + '/*.' + yyyymmdd + '.nc4')

        # Download data 
        print('Download ' + yyyymmdd)
        print('Directory is ' + data_dir)
        get_tavg1(yyyymmdd, data_dir, collections=collections,
                user=user, password=password)

        # Next day
        nextDate_D = currDate_D + datetime.timedelta(days=1)
        if (str(nextDate_D)[5:7] != mm):
            break
        currDate_D = nextDate_D
#
#------------------------------------------------------------------------------
#
def get_tavg1_time_range(startDate, endDate, root_data_dir,
        collections=None,
        user='yi-wang-4@uiowa.edu',
        password=''):
    """ Download GEOS-FP asm tavg1 (collecations is None)data 
    of a time range.
    But users can specific variables now.
    Parameters
    ----------
    startDate : str
        Start date, '20180501' for example
    endDate : str
        End date, '20180531' for example
    root_data_dir : str
        The root directory where all data will be saved to.
    collections : list
        list of variables.
    user : str
        MERRA-2 website user name.
    password : str
        MERRA-2 website password.

    """

    if root_data_dir[-1] != '/':
        root_data_dir = root_data_dir + '/'

    # Date
    currDate   = startDate
    currDate_D = datetime.datetime.strptime(currDate, '%Y%m%d')
    endDate_D  = datetime.datetime.strptime(endDate,  '%Y%m%d')

    # loop everyday
    while currDate_D <= endDate_D:
        
        # Date
        currDate = str(currDate_D)
        yyyy = currDate[0:4]
        mm   = currDate[5:7]
        dd   = currDate[8:10]
        yyyymmdd = yyyy + mm + dd

        # Directory
        data_dir = root_data_dir + yyyy + '/' + mm
        if not os.path.isdir(data_dir):
            os.system('mkdir -p ' + data_dir)
        os.system('rm -f ' + data_dir + '/*.' + yyyymmdd + '.nc4')

        # Download data 
        print('Download ' + yyyymmdd)
        print('Directory is ' + data_dir)
        get_tavg1(yyyymmdd, data_dir, collections=collections,
                user=user, password=password)

        # go to next day
        currDate_D = currDate_D + datetime.timedelta(days=1)
#
#------------------------------------------------------------------------------
#
