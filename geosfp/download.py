"""
Created on March 17, 2020

@author: Yi Wang
"""

import datetime
import os


#
#------------------------------------------------------------------------------
#
def tavg1_url(yyyymmdd):
    """ Generate a list that contains the URLs of GEOS-FP 
    one-day asm tavg1 data.
    (ywang, 03/17/2020)

    Parameters
    ----------
    yyyymmdd : str
        Date, '20180501' for example

    Returns
    -------
    urls : list
        All URLs
    """

    # URL of GEOS-FP data
    url = 'https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/'

    #
    collections = ['tavg1_2d_flx_Nx', 'tavg1_2d_lnd_Nx',
                   'tavg1_2d_rad_Nx', 'tavg1_2d_slv_Nx']

    # Generate URLs
    urls = []
    yyyy = yyyymmdd[0:4]
    mm   = yyyymmdd[4:6]
    dd   = yyyymmdd[6:8]
    url_curr_day = url + 'Y' + yyyy + '/M' + mm + '/D' + dd + '/'
    Nh = 24
    for ih in range(Nh):

        ch = str(ih).zfill(2) + '30'

        for coln in collections:

            filename = url + 'Y' + yyyy + '/M' + mm + '/D' + dd + \
                    '/GEOS.fp.asm.' + coln  + '.' + yyyymmdd  + '_' + \
                    ch + '.V01.nc4'

            urls.append(filename)

    return urls
#
#------------------------------------------------------------------------------
#
def get_tavg1(yyyymmdd, data_dir, retry=10):
    """ Download GEOS-FP one-day asm tavg1 data.
    (ywang, 03/17/2020)

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

    Returns
    -------
    out_dict : dict

    """

    out_dict = {}

    # Save filename to be downloaded to tmp_tavg1_files.txt
    tmp_file = './tmp_tavg1_files.txt'
    f = open(tmp_file, 'w')
    urls = tavg1_url(yyyymmdd)
    out_dict['urls'] = urls
    for filename in urls:
        f.write(filename + '\n')

    # Close tmp_tavg1_files.txt
    f.close()

    # Download data
    cmd = 'wget -i ' + tmp_file + ' -P ' + data_dir
    os.system(cmd)

    # Remove tmp_tavg1_files.txt
    os.system('rm -rf ' + tmp_file)

    # check
    check_tavg1(urls, data_dir, retry=retry)

    return out_dict
#
#------------------------------------------------------------------------------
#
def check_tavg1(yyyymmdd, data_dir, retry=10):
    """ Check if any GEOS-FP one-day asm tavg1 data are not
    downloaded successfully. If so, download them again.
    (ywang, 03/17/2020)

    Parameters
    ----------
    yyyymmdd : str or list
        If str, Date, '20180501' for example
        If list, it contains the URLs of GEOS-FP 
        one-day asm tavg1 data.
    data_dir : str
        The directory where all data will be saved to.

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
        urls = tavg1_url(yyyymmdd)
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
                wget_cmd = 'wget ' + url + ' -O ' + filename
                print('Retry "' + wget_cmd + '" {:} time(s).'.format(itry))
                os.system(wget_cmd)

        # Fail after maxium retry.
        if (os.system(nc_cmd) != 0) and (itry == retry):
            print(url + ' failed after {} tries'.format(itry))
#
#------------------------------------------------------------------------------
#
def get_tavg1_month(yyyymm, root_data_dir):
    """ Download GEOS-FP one-month asm tavg1 data.
    (ywang, 03/18/2020)

    Parameters
    ----------
    yyyymm : str
        Month, '201805' for example.
    root_data_dir : str
        The root directory where all data will be saved to.

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
        data_dir = root_data_dir + 'Y' + yyyy + '/M' + mm + '/D' + dd
        if not os.path.isdir(data_dir):
            os.system('mkdir -p ' + data_dir)
        os.system('rm -f ' + data_dir + '/*.nc4')

        # Download data 
        print('Download ' + yyyymmdd)
        print('Directory is ' + data_dir)
        get_tavg1(yyyymmdd, data_dir)

        # Next day
        nextDate_D = currDate_D + datetime.timedelta(days=1)
        if (str(nextDate_D)[5:7] != mm):
            break
        currDate_D = nextDate_D
#
#------------------------------------------------------------------------------
#
def get_tavg1_time_range(startDate, endDate, root_data_dir):
    """ Download GEOS-FP asm tavg1 data of a time range.
    Parameters
    ----------
    startDate : str
        Start date, '20180501' for example
    endDate : str
        End date, '20180531' for example
    root_data_dir : str
        The root directory where all data will be saved to.
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
        data_dir = root_data_dir + 'Y' + yyyy + '/M' + mm + '/D' + dd
        if not os.path.isdir(data_dir):
            os.system('mkdir -p ' + data_dir)
        os.system('rm -f ' + data_dir + '/*.nc4')

        # Download data 
        print('Download ' + yyyymmdd)
        print('Directory is ' + data_dir)
        get_tavg1(yyyymmdd, data_dir)

        # go to next day
        currDate_D = currDate_D + datetime.timedelta(days=1)
#
#------------------------------------------------------------------------------
#
