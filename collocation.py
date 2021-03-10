"""
Created on September 4, 2019

@author: Yi Wang
"""

import datetime
from geopy.distance import great_circle
import numpy as np
import pandas as pd

from mylib.select_region import is_in_the_granule


def find_aeronet_granule(sat_lat, sat_lon, aeronet_site):
    """ Find all AERONET sites that are in the rectangle that
    encompass the granule.

    Parameters
    ----------
    sat_lat :
        Satellite latitude
    sat_lon :
        Satellite longitude
    aeronet_site : 
        Data returned by read_aeronet_site()

    Returns
    -------
    out_data : pd.DataFrame
        All sites found are store in the out_data.
        Keys are 'site_name', 'site_lat', and 'site+lon'.

    """

    # AERONET site names, latitudes, and longitudes
    site_name = np.array(aeronet_site['Site'])
    site_lat  = np.array(aeronet_site['Latitude'])
    site_lon  = np.array(aeronet_site['Longitude'])

    # Find all site
    flag = is_in_the_granule(sat_lat, sat_lon, site_lat, site_lon)

    # save all site information to a dict
    out_data = {}
    if np.any(flag):
        out_data['site_name'] = site_name[flag]
        out_data['site_lat']  = site_lat[flag]
        out_data['site_lon']  = site_lon[flag]

    out_data = pd.DataFrame(out_data)

    return out_data

def collocate_site_satellite(site_TAI93, site_var, site_lat, site_lon,
        sat_TAI93, sat_lat, sat_lon, sat_var_dict, diameter=50.0,
        window=None, sat_var_min=-0.05):
    """ Calculate satellite overpass time, then calculate aeronet
    measurement average at the overpass time using data within 
    a time window and the average of satellite variables within
    a circle with diameter of *diameter* and centered at site

    Parameters
    ----------
    site_TAI93 : 1-D array, float
        Site observational time. Seconds since 1993-01-01T00:00:00
    site_var : 1-D array, float
        Site measurements. It has same dimension as site_TAI93
    site_lat : float
        Site latitude
    site_lon : float
        Site longitude
    sat_TAI93 : array, float
        Satellite observational time. Seconds since 1993-01-01T00:00:00
    sat_lat : array, float
        Satellite latitude
    sat_lon : array, float
        Satellite longitude
    sat_var_dict : dict
        The vaules are satellite variables.
    diameter : float
        Diameter of the circle that is centered at site.
        Unit is km.
    window : float
        Time window used for calculate average. Unit is second.
        If its vaule is None, the vaule is assigned to
        (diameter / 50.0) * 3600.0
    sat_res : float
        Satellite spatial resolution (km) at nadir. If the along track
        and cross track spatical resolution are unequal, just assign
        the maximum. This parameter is used to save compuational
        time. If None is assigned, the distance between every pixel
        and site will be calcualted.
    sat_var_min : float
        Satellite variable that is less then *sat_var_min* is
        undefined.

    Returns
    -------
    TAI93 : float or None
        Average of satellite overpass time
        If TAI93 is None, it means the satellie didn't overpass
        the site
    (site_ave, site_N) : (float, float)
        Average of *site_var* and number of aeronet measurements
        within the time window.
    sat_result : dict
        Keys are the keys from sat_var_dict. Vaules are tuples of
        average of satellite variable and number of variables
        used to calculate average.

    """

    # Initialize return values
    TAI93 = None
    site_ave = None
    site_N = None
    sat_result = {}
    for sat_var_name in sat_var_dict:
        sat_result[sat_var_name] = (None, 0)

    # Time winodw
    # Average wind speed is 50 km/h or 50 km per 3600 seconds.
    if window is None:
        window = (diameter / 50.0) * 3600.0

    half_window = window / 2.0

    # Radius of the circle that is centered as AERONET site.
    radius = diameter / 2.0

    # Shrink searching areas if possible.
    sat_TAI93_region = np.array(sat_TAI93)
    sat_lat_region   = np.array(sat_lat)
    sat_lon_region   = np.array(sat_lon)

    deg_distance_thre = radius / 10.0

    min_lat = site_lat - deg_distance_thre
    max_lat = site_lat + deg_distance_thre
    min_lon = site_lon - deg_distance_thre
    max_lon = site_lon + deg_distance_thre

    flag_lat = np.logical_and(sat_lat > min_lat, sat_lat < max_lat)
    flag_lon = np.logical_and(sat_lon > min_lon, sat_lon < max_lon)
    flag_region = np.logical_and(flag_lat, flag_lon)

    # Process data in the possible region
    if np.any(flag_region):

        sat_TAI93_region = sat_TAI93_region[flag_region]
        sat_lat_region   = sat_lat_region[flag_region]
        sat_lon_region   = sat_lon_region[flag_region]


        # Find satellite data in the circle
        site_lat_lon = (site_lat, site_lon)
        flag_circle = np.empty(sat_TAI93_region.shape, dtype=bool)
        for i in range(len(sat_TAI93_region)):

            # latitude and longitude of a satellite pixel
            sat_pixel_lat_lon = (sat_lat_region[i], sat_lon_region[i])

            # great circle distance
            gc_discance = great_circle(site_lat_lon, sat_pixel_lat_lon)

            # Pixel in the circle? 
            if (gc_discance <= radius):
                flag_circle[i] = True
            else:
                flag_circle[i] = False

        # Process data in the circle 
        if np.any(flag_circle):

            # get all satellite TAI93, lat, and lon in the circle
            sat_TAI93_circle = sat_TAI93_region[flag_circle]
            sat_lat_circle   = sat_lat_region[flag_circle]
            sat_lon_circle   = sat_lon_region[flag_circle]

            # calculate satellite overpass time
            TAI93 = np.mean(sat_TAI93_circle)

            # calculate time window
            TAI93_beg = TAI93 - half_window
            TAI93_end = TAI93 + half_window
            flag_window = np.logical_and((site_TAI93 > TAI93_beg),
                                         (site_TAI93 < TAI93_end))

            # Number of site observations within time window
            site_N = np.sum(flag_window)
            if site_N > 0:

                # calculate site mean
                site_ave = np.mean(site_var[flag_window])

                # get all satellite variables in the circle
                sat_var_circle_dict = {}
                for sat_var_name in sat_var_dict:
                    sat_var_circle_dict[sat_var_name] = \
                        sat_var_dict[sat_var_name][flag_region][flag_circle]

                # calculate satellite average in the circle
                for sat_var_name in sat_var_circle_dict:

                    # satellite variable in the circle
                    var = sat_var_circle_dict[sat_var_name]

                    # flag for meaningful vaules
                    flag_var = (var >= sat_var_min)

                    var_N = np.sum(flag_var)
                    if var_N > 0:

                        # satellite average in the circle
                        var_ave = np.mean(var[flag_var])

                        sat_result[sat_var_name] = (var_ave, var_N)

                    else:

                        sat_result[sat_var_name] = (-999.0, var_N)

    return TAI93, (site_ave, site_N), sat_result

def collocate_line_satellite(ori_line_TAI93, ori_line_var, 
        ori_line_lat, ori_line_lon,
        sat_TAI93, sat_lat, sat_lon, sat_var_dict,
        collocation_result, name,
        obs_var_name='man_aod550', obs_var_N_name='man_aod550_N',
        diameter=50.0, window=None, sat_var_min=-0.05,
        ):
    """ Collocate "line" (For example, ship) observations 
    with satellite retrievals.
    (Yi Wang, 12/01/2019)

    Parameters
    ----------
    ori_line_TAI93 : 1-D array, float
        "Line" observational time. Seconds since 1993-01-01T00:00:00
    ori_line_var : 1-D array, float
        "Line" measurements. It has same dimension as ori_line_TAI93
    ori_line_lat : float
        "Line" latitude
    ori_line_lon : float
        "Line" longitude
    sat_TAI93 : array, float
        Satellite observational time. Seconds since 1993-01-01T00:00:00
    sat_lat : array, float
        Satellite latitude
    sat_lon : array, float
        Satellite longitude
    sat_var_dict : dict
        The vaules are satellite variables.
    collocation_result : dict
        collocation_result parameter for save_collocation_to_dict
    name : str
        site_name parameter for save_collocation_to_dict
    obs_var_name : str
        site_var_name parameter for save_collocation_to_dict
    obs_var_N_name : str
        site_var_N_name parameter for save_collocation_to_dict
    diameter : float
        Diameter of the circle that is centered at site.
        Unit is km.
    window : float
        Time window used for calculate average. Unit is second.
        If its vaule is None, the vaule is assigned to
        (diameter / 50.0) * 3600.0
    sat_var_min : float
        Satellite variable that is less then *sat_var_min* is
        undefined.

    Returns
    -------
    out_dict : dict
        key1: 'collocated' (logical)
        key2: 'collocation_result' (dict), if collocated is T.
        *collocation_result* is also modified

    """

    out_dict = {}

    # Time window
    # Average wind speed is 50 km/h or 50 km per 3600 seconds.
    if window is None:
        window = (diameter / 50.0) * 3600.0

    # half time window
    half_window = window / 2.0

    # As we can only determine the range of of referece TAI93,
    # we rougthly estimate TAI_beg and TAI_end
    sat_sec_range = np.max(sat_TAI93) - np.min(sat_TAI93)
    sat_sec_mean = np.mean(sat_TAI93)
    TAI_beg = sat_sec_mean - half_window - sat_sec_range
    TAI_end = sat_sec_mean + half_window + sat_sec_range

    # get a subset accroding to TAI93 
    flag = np.logical_and(ori_line_TAI93>=TAI_beg, ori_line_TAI93<=TAI_end)
    if np.sum(flag) > 0:
        line_TAI93 = ori_line_TAI93[flag]
        line_var   = ori_line_var[flag]
        line_lat   = ori_line_lat[flag]
        line_lon   = ori_line_lon[flag]
    else:
        out_dict['collocated'] = False
        return out_dict

    # collocate and save
    for i in range(len(line_TAI93)):
        one_line_TAI93 = np.array([line_TAI93[i]])
        one_line_var   = np.array([line_var[i]])
        one_line_lat   = line_lat[i]
        one_line_lon   = line_lon[i]
        # collocate
        TAI93, (site_ave, site_N), sat_result = \
                collocate_site_satellite(one_line_TAI93, one_line_var,
                        one_line_lat, one_line_lon,
                        sat_TAI93, sat_lat, sat_lon, sat_var_dict,
                        diameter=diameter, window=window,
                        sat_var_min=sat_var_min)

        # save data to dictionary 
        sat_flag = False
        for sat_key in sat_result:
            if sat_result[sat_key][1] > 0:
                # at least one satellite variable
                sat_flag = True
        if (TAI93 is not None) and (site_ave is not None) and sat_flag:

            save_collocation_to_dict(
                    collocation_result, TAI93, site_ave, site_N,
                    sat_result, one_line_lat, one_line_lon, name,
                    site_var_name=obs_var_name,
                    site_var_N_name=obs_var_N_name
                    )

            out_dict['collocated'] = True
            out_dict['collocation_result'] = collocation_result

    return out_dict




def save_collocation_to_dict(collocation_result, TAI93, site_var, site_var_N,
        sat_result, lat, lon, site_name, site_var_name='site_aod550', 
        site_var_N_name='site_aod550_N'):
    """ Save collocation to dictionary

    Parameters
    ----------
    collocation_result : dict
        collocation_result[site_name] = {
                'ymd'   : [],
                'hms'   : [],
                'TAI93' : [],
                site_var_name  : [],
                site_var_N_name : [],
                sat_result.keys() : [],
                sat_result.keys()+'_N' : [],
                'lat' : []
                'lon' : []
                }

    Returns
    -------
    None

    """

    # convert TAI93 to YYYY-MM-DD HH:MM:SS
    ref_time = datetime.datetime.strptime('1993-01-01 00:00:00',
                                          '%Y-%m-%d %H:%M:%S')
    curr_time = ref_time + datetime.timedelta(seconds=round(TAI93))
    ymd = str(curr_time).split()[0]
    hms = str(curr_time).split()[1]

    # ymd, hms, and TAI93
    collocation_result[site_name]['ymd'].append(ymd)
    collocation_result[site_name]['hms'].append(hms)
    collocation_result[site_name]['TAI93'].append(TAI93)

    # site_var, site_var_N
    collocation_result[site_name][site_var_name].append(site_var)
    collocation_result[site_name][site_var_N_name].append(site_var_N)

    # sat_result
    #print(sat_result)
    #print(site_name)
    #print(collocation_result[site_name].keys())
    for key in sat_result:
        sat_var   = sat_result[key][0]
        sat_var_N = sat_result[key][1]
        collocation_result[site_name][key].append(sat_var)
        collocation_result[site_name][key+'_N'].append(sat_var_N)

    # lat and lon
    collocation_result[site_name]['lat'].append(lat)
    collocation_result[site_name]['lon'].append(lon)

    return None
