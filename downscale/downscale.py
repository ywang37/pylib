"""
Created on Feburary 15, 2020

@author: Yi Wang
"""

import numpy as np
import sys
from tqdm import tqdm

from mylib.downscale.resample_parallel import alber_equal_area
from mylib.downscale.resample_parallel import getWeights

#
#------------------------------------------------------------------------------
#
def search_possible_overlap_grid(limit_src, lat_e_tar, lon_e_tar):
    """ Find all target grids that overlap with source limit region

    Parameters
    ----------
    limit_src : 1-D array-like
        Source region limit
        [lat_min, lon_min, lat_max, lon_max]
    lat_e_tar: 2-D array
        Latitude edges of 2-D target grid
    lon_e_tar: 2-D array
        Longitude edges of 2-D target grid

    Returns
    ind_dict :
        Keys:
            find: bool variabe. True and False means find
                and not find possible overlap grid, respectively.
            lat_i0: target grid latitude start index
            lat_i1: target grid latitude end index
            lon_j0: target grid longitude start index
            lon_j1: target grdi longitude end index
    """

    src_lat_min = limit_src[0]
    src_lat_max = limit_src[2]
    src_lon_min = limit_src[1]
    src_lon_max = limit_src[3]

    # possible overlap grid edges
    flag_lat = np.logical_and(lat_e_tar >= src_lat_min, 
            lat_e_tar <= src_lat_max)
    flag_lon = np.logical_and(lon_e_tar >= src_lon_min, 
            lon_e_tar <= src_lon_max)
    flag = np.logical_and(flag_lat, flag_lon)

    ind_dict = {}
    if np.sum(flag) == 0:
        ind_dict['find'] = False
        return ind_dict
    else:
        ind_dict['find'] = True

    # target dimension
    dim_tar = lat_e_tar.shape
    nlate = dim_tar[0]
    nlone = dim_tar[1]

    # find lat_i0: target grid latitude start index
    for i0 in range(nlate):
        if np.any(flag[i0,:]):
            ind_dict['lat_i0'] = i0
            if ind_dict['lat_i0'] > 0:
                ind_dict['lat_i0'] -= 1
            break

    # find lat_i1: target grid latitude end index
    for i1 in np.array(range(nlate))[::-1]:
        if np.any(flag[i1,:]):
            ind_dict['lat_i1'] = i1
            if ind_dict['lat_i1'] < (nlate-1):
                ind_dict['lat_i1'] += 1
            break

    # find lon_j0: target grid longitude start index
    for j0 in range(nlone):
        if np.any(flag[:,j0]):
            ind_dict['lon_j0'] = j0
            if ind_dict['lon_j0'] > 0:
                ind_dict['lon_j0'] -= 1
            break

    # find lon_j1: target grdi longitude end index
    for j1 in np.array(range(nlone))[::-1]:
        if np.any(flag[:,j1]):
            ind_dict['lon_j1'] = j1
            if ind_dict['lon_j1'] < (nlone-1):
                ind_dict['lon_j1'] += 1
            break

    return ind_dict
#
#------------------------------------------------------------------------------
#
def calc_overlap_area_record(lat_e_src, lon_e_src, lat_e_tar, lon_e_tar,
        flag_src=None):
    """ Calculate overlapping area of source grids and
    target grids

    Parameters
    ----------
    lat_e_src: 2-D array
        Latitude edges of 2-D source grid
    lon_e_src: 2-D array
        Longitude edges of 2-D source grid
    lat_e_tar: 2-D array
        Latitude edges of 2-D target grid
    lon_e_tar: 2-D array
        Longitude edges of 2-D target grid
    flag_src: 2-D array
        Dimension length is 1 smaller than that of *lat_e_src*
        Only process grid with True flag.
        If flag_src is None, all grids are processed.

    Returns
    -------
    src_overlap_area : dict
        Overlapping area
    tar_overlap_ind : dict
        index

    """

    # source dimension
    dim_e_src = lat_e_src.shape
    dim_c_src = (dim_e_src[0] - 1, dim_e_src[1] - 1)

    # target dimension
    dim_e_tar = lat_e_tar.shape
    dim_c_tar = (dim_e_tar[0] - 1, dim_e_tar[1] - 1)

    # source flag
    if flag_src is None:
        flag_src = np.full(dim_c_src, True)
    else:
        if (flag_src.shape != dim_c_src):
            print(' - calc_overlap_area_record: flag_src dimension error.')
            print('flag_src dimension is: ', flag_src.shape)
            print('dim_c_src dimension is: ', dim_c_src)
            exit()

    # Convert to alber equal area
    lon_0 = 110.0
    x_e_src, y_e_src = alber_equal_area(lon_e_src, lat_e_src, lon_0=lon_0)
    x_e_tar, y_e_tar = alber_equal_area(lon_e_tar, lat_e_tar, lon_0=lon_0)

    # dict to stroe overlap area
    src_overlap_area = {}

    # dict to store overlap source grid indexes for target grids
    tar_overlap_ind = {}
    for i_tar in range(dim_c_tar[0]):
        for j_tar in range(dim_c_tar[1]):
            tar_overlap_ind[(i_tar,j_tar)] = []

    #
    print(' - calc_overlap_area_record()')
    for i_src in tqdm(range(dim_c_src[0])):
        for j_src in range(dim_c_src[1]):

            # only process data with flag is True
            if  (not flag_src[i_src,j_src]):
                continue

            # find all target grids that overlap with source grid
            lat_corner = [lat_e_src[i_src,j_src], lat_e_src[i_src,j_src+1],
                    lat_e_src[i_src+1,j_src+1], lat_e_src[i_src+1,j_src]]
            lon_corner = [lon_e_src[i_src,j_src], lon_e_src[i_src,j_src+1],
                    lon_e_src[i_src+1,j_src+1], lon_e_src[i_src+1,j_src]]
            limit_src = [min(lat_corner), min(lon_corner),
                    max(lat_corner), max(lon_corner)]
            tar_ind_dict = search_possible_overlap_grid(limit_src, 
                    lat_e_tar, lon_e_tar)

            # no overlap
            if not tar_ind_dict['find']:
                continue
            #print(limit_src)
            i0_tar = tar_ind_dict['lat_i0']
            i1_tar = tar_ind_dict['lat_i1']
            j0_tar = tar_ind_dict['lon_j0']
            j1_tar = tar_ind_dict['lon_j1']
            #print(lat_e_tar[i0_tar:i1_tar+1,j0_tar:j1_tar+1])
            #print(lon_e_tar[i0_tar:i1_tar+1,j0_tar:j1_tar+1])
            #print(tar_ind_dict)
            #print('i0_tar = {}'.format(i0_tar))
            #print('i1_tar = {}'.format(i1_tar))
            #print('j0_tar = {}'.format(j0_tar))
            #print('j1_tar = {}'.format(j1_tar))

            # cPixel for calculate overlapping area
            cPixel = {}
            cPixel['Upper_Left_Lon']  = x_e_src[i_src+1,j_src]
            cPixel['Upper_Right_Lon'] = x_e_src[i_src+1,j_src+1]
            cPixel['Lower_Right_Lon'] = x_e_src[i_src,j_src+1]
            cPixel['Lower_Left_Lon']  = x_e_src[i_src,j_src]
            cPixel['Upper_Left_Lat']  = y_e_src[i_src+1,j_src]
            cPixel['Upper_Right_Lat'] = y_e_src[i_src+1,j_src+1]
            cPixel['Lower_Right_Lat'] = y_e_src[i_src,j_src+1]
            cPixel['Lower_Left_Lat']  = y_e_src[i_src,j_src]

            #print(cPixel)

            # surPixels for calculate overlapping area
            surPixels = {}
            surPixels['Upper_Left_Lon'] = \
                    x_e_tar[i0_tar+1:i1_tar+1,j0_tar:j1_tar]
            surPixels['Upper_Right_Lon'] = \
                    x_e_tar[i0_tar+1:i1_tar+1,j0_tar+1:j1_tar+1]
            surPixels['Lower_Right_Lon'] = \
                    x_e_tar[i0_tar:i1_tar,j0_tar+1:j1_tar+1]
            surPixels['Lower_Left_Lon'] = \
                    x_e_tar[i0_tar:i1_tar,j0_tar:j1_tar]
            surPixels['Upper_Left_Lat'] = \
                    y_e_tar[i0_tar+1:i1_tar+1,j0_tar:j1_tar]
            surPixels['Upper_Right_Lat'] = \
                    y_e_tar[i0_tar+1:i1_tar+1,j0_tar+1:j1_tar+1]
            surPixels['Lower_Right_Lat'] = \
                    y_e_tar[i0_tar:i1_tar,j0_tar+1:j1_tar+1]
            surPixels['Lower_Left_Lat'] = \
                    y_e_tar[i0_tar:i1_tar,j0_tar:j1_tar]

            #print(surPixels)
            

            # calculate overlapping area
            weights, area, cArea = getWeights(cPixel, surPixels, None)
#            print('weight')
#            print(weights)
#            print('area')
#            print(area)
#            print('cArea')
#            print(cArea)
#            print('sum of area')
#            print(np.sum(area))
#            print('sum of weight')
#            print(np.sum(weights))
#            exit()

            # dict to store overlap area
            src_overlap_area[(i_src,j_src)] = {}
            src_overlap_area[(i_src,j_src)]['area'] = cArea
            for i in range(area.shape[0]):
                for j in range(area.shape[1]):

                    # target index
                    i_tar = i + i0_tar
                    j_tar = j + j0_tar

                    # 
                    src_overlap_area[(i_src,j_src)][(i_tar,j_tar)] = \
                            area[i,j]

                    tar_overlap_ind[(i_tar,j_tar)].append((i_src,j_src))

            #print(src_overlap_area[(i_src,j_src)])


            #for i_tar in range(dim_c_tar[0]):
            #    for j_tar in range(dim_c_tar[1]):
            #        if len(tar_overlap_ind[(i_tar,j_tar)]) > 0:
            #            print((i_tar,j_tar), tar_overlap_ind[(i_tar,j_tar)])

            #exit()


    return src_overlap_area, tar_overlap_ind
#
#------------------------------------------------------------------------------
#
def downscale_emissions(src_overlap_area, tar_overlap_ind, src_emi, tar_emi,
        src_format):
    """ Downscale source emissions *src_emi* according to the
    distribution of target emisisons *tar_emi*

    Parameters
    ----------
    src_overlap_area : dict
        Overlapping area from function calc_overlap_area_record
    tar_overlap_ind : dict
        index from calc_overlap_area_record
    src_emi : 2-D array
        Source emissions
    tar_emi : 2-D array
        Target emissions (It must be emissions rates; how much
        emissions in unit area)
    src_format : str
        'total' : emisisons of every sources grid
        'rate'  : emisisons rates of every sources grid

    Returns
    -------
    out_emi
    out_overlap_area

    """

    # downscaled emissions at target grids
    out_emi = np.zeros_like(tar_emi)

    # overlap area at target grids
    out_overlap_area = np.zeros_like(tar_emi)

    # distribute source emissions *src_emi* according to target grid
    src_emi_dist = {}
    for src_ij in src_overlap_area:

        # one source grid
        src = src_overlap_area[src_ij]

        # calculate weight
        weight = {}
        tmp_sum = 0.0
        for tar_ij in src:
            if tar_ij != 'area':
                weight[tar_ij] = src[tar_ij] * tar_emi[tar_ij]
                tmp_sum += weight[tar_ij]
        if tmp_sum > 0.0:
            for tar_ij in weight:
                weight[tar_ij] /= tmp_sum
        else:
            for tar_ij in weight:
                weight[tar_ij] = src[tar_ij] / src['area']
        for tar_ij in weight:
            if weight[tar_ij] < 0.0:
                print(' - downscale_emissions: weight error.')
                exit()

        # distribute source emissions
        src_emi_dist[src_ij] = {}
        if src_format == 'rate':
            one_src_emi = src_emi[src_ij] * src['area']
        elif src_format == 'total':
            one_src_emi = src_emi[src_ij]
        else:
            print(' - downscale_emissions: src_format error.')
            print('src_format: ' + src_format)
            exit()
        for tar_ij in weight:
            src_emi_dist[src_ij][tar_ij] = one_src_emi * weight[tar_ij]

    # aggregate source emissions to target grid
    for i_tar in range(out_emi.shape[0]):
        for j_tar in range(out_emi.shape[1]):

            # get index list
            one_tar_overlap_ind = tar_overlap_ind[(i_tar, j_tar)]

            # aggregate emissions
            if len(one_tar_overlap_ind) > 0:
                for i_o in range(len(one_tar_overlap_ind)):
                    # source grid index
                    src_ij = one_tar_overlap_ind[i_o]
                    # aggregate emissions
                    out_emi[i_tar,j_tar] += \
                            src_emi_dist[src_ij][(i_tar,j_tar)]
                    # aggregate area
                    out_overlap_area[i_tar,j_tar] += \
                            src_overlap_area[src_ij][(i_tar,j_tar)]

    # change back to rate
    if src_format == 'rate':
        area_flag = np.logical_and( out_overlap_area > 0.0, out_emi > 0.0)
        out_emi[area_flag] = out_emi[area_flag] / out_overlap_area[area_flag]

    return out_emi, out_overlap_area




                    
                    














#
#------------------------------------------------------------------------------
#
