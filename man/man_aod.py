"""
Created on November 25, 2019

@author: Yi Wang
"""

from astropy.time import Time
import glob
import numpy as np
import pandas as pd
import sys

from mylib.aeronetpy import readAeronet
from mylib.pro_aerosol import cal_ae, interpolate_aod

sys.path.append('/Dedicated/jwang-data/ywang/project_global_CWA/shared_code')
from aeronet_io import read_aeronet_site

#
#------------------------------------------------------------------------------
#
def get_man_550_aod(in_file, out_dir=None):
    """ Get MAN 550 nm AOD.
    """

    varnames = [
            'Latitude',
            'Longitude',
            'AOD_440nm',
            'AOD_500nm',
            'AOD_675nm',
            ]

    # output file
    filename = in_file.split('/')[-1]
    filename = filename.replace('.', '_')
    filename = filename + '.csv'
    if out_dir is not None:
        if out_dir[-1] != '/':
            out_dir = out_dir + '/'
        out_file = out_dir + filename

    #
    if (('_daily_' in filename) or ('_series_' in filename)):
        nobs_flag = True
    else:
        nobs_flag = False

    print(' - get_man_550_aod: Process ' + in_file)

    # read aeronet data
    if nobs_flag:
        varnames.append('Number_of_Observations')
    data = readAeronet(in_file, vars=varnames, na_values=None)

    # no data
    if isinstance(data, list):
        print(" - get_man_550_aod WARNING: There is no data.")
        return None

    # replace nan by -999.0
    data.pop('DailyIndex')
    data = pd.DataFrame(data)
    data = data.fillna(-999.0)

    ymd = np.array(data['Date'],dtype=str)
    hms = np.array(data['Time'],dtype=str)
    aod440 = data['AOD_440nm']
    aod500 = data['AOD_500nm']
    aod675 = data['AOD_675nm']
    lat = data['Latitude']
    lon = data['Longitude']
    if nobs_flag:
        nobs = data['Number_of_Observations']

    # filter data
    flag440 = np.logical_and(aod440>0, aod675>0)
    flag500 = np.logical_and(aod500>0, aod675>0)
    if ( (np.sum(flag440) == 0) and (np.sum(flag500) == 0) ):
        print(" - get_man_550_aod WARNING: " + \
                "There is no data for interpolation.")
        return None
    if ( np.sum(flag440) >= np.sum(flag500) ):
        flag = flag440
        ref_band = 440
    else:
        flag = flag500
        ref_band = 500

    ymd = ymd[flag]
    hms = hms[flag]
    aod440 = aod440[flag]
    aod500 = aod500[flag]
    aod675 = aod675[flag]
    lat = lat[flag]
    lon = lon[flag]
    if nobs_flag:
        nobs = nobs[flag]

    if ref_band == 440:

        # calculate angstrom exponent
        AE_calculated = cal_ae(440.0, 675.0, aod440, aod675)
        AE_band = 'AE_440-675'

        # interpolate aod550
        aod550 = interpolate_aod(440.0, aod440, AE_calculated, 550.0)

    elif ref_band == 500:

        # calculate angstrom exponent
        AE_calculated = cal_ae(500.0, 675.0, aod500, aod675)
        AE_band = 'AE_500-675'

        # interpolate aod550
        aod550 = interpolate_aod(500.0, aod500, AE_calculated, 550.0)

    else:
        print(" - get_man_550_aod: !!! ERROR !!!")
        exit()

    # convert date to seconds since 1993-1-1 00:00:00
    ymdhms = np.core.defchararray.add(ymd,
                np.core.defchararray.add('T', hms))
    one_day_TAI = 24 * 3600E0
    TAI93 = (Time(ymdhms).jd - Time('1993-01-01T00:00:00').jd) * one_day_TAI

    # save to csv
    df = {
        'ymd': ymd,
        'hms': hms,
        'TAI93': TAI93,
        'aod440': aod440,
        'aod500': aod500,
        'aod550': aod550,
        'aod675': aod675,
#        AE_band: AE_calculated,
        'lat': lat,
        'lon': lon,
        }
    df = pd.DataFrame(df)
    if nobs_flag:
        df['nobs'] = nobs

    if out_dir is not None:
        print(' - Output ' + out_file)
        df.to_csv(out_file, index=False)

    return df
#
#------------------------------------------------------------------------------
#
