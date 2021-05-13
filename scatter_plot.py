import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import mylib.metrics as ms

def histogram2D(ax, in_x_data, in_y_data, bins=10, range=None,
        weights=None, density=None, 
        cmap=plt.get_cmap('jet'), vmin=0, vmax=None, alpha=0.7,
        **kwargs):
    """ plot 2D histogram of density

    Parameters:
    bins, range, weights, and density are passed to
        np.histogram2d
    cmap :
        colormap passed ax.pcolormesh
    vmin : 
        The minumum valule that is colored
    vmax :
        The maximum value that is colored
    alpha :
        Tranparency
    **kwargs is passed to ax.pcolormesh

    returns
    -------
    mesh : The object returned by ax.pcolormesh

    """

    # copy and flatten data
    x_data = in_x_data.flatten()
    y_data = in_y_data.flatten()

    # get a mask with those elements posible to compare (non-nan values)
    mask = np.logical_and(np.logical_not(np.isnan(x_data)), 
                          np.logical_not(np.isnan(y_data)))
    # useful data
    x_data = x_data[mask]
    y_data = y_data[mask]

    # statistics of 2D histogram
    H, xedges, yedges = np.histogram2d(x_data, y_data,
            bins=bins, range=range, weights=weights, density=density)
    H = H.T # Let each row list bins with common y range

    H = np.ma.masked_array(H, H<vmin)


    # X and Y coordinate
    X, Y = np.meshgrid(xedges, yedges)

    mesh = ax.pcolormesh(X, Y, H, cmap=cmap, 
            vmin=vmin, vmax=vmax, alpha=alpha)

    return mesh

def scatter(ax, in_x_data, in_y_data, \
        label_ul=None, label_lr=None,
        text_color=None,
        regress_line=True, one2one_line=True, **kwargs):
    """ Function to create scatter plots.
    """

    if label_ul is None:
        label_ul = [
                'R',
                'x_mean_std',
                'y_mean_std',
                'rmse',
                'N'
                ]
    
    # copy and flatten data
    x_data = in_x_data.flatten()
    y_data = in_y_data.flatten()

    # get a mask with those elements posible to compare (non-nan values)
    mask = np.logical_and(np.logical_not(np.isnan(x_data)), 
                          np.logical_not(np.isnan(y_data)))
    n_colocations = len(mask[mask==True])
    x_data = x_data[mask]
    y_data = y_data[mask]
    if n_colocations == 0:
        return
    
    # liner regression
    slope, intercept, correlation, p_value_slope, std_error = \
            stats.linregress(x_data, y_data)

    # Calculates a Pearson correlation coefficient 
    # and the p-value for testing non-correlation
    r, p_value = stats.pearsonr(x_data, y_data)
    r_s = "R = {:.2f}".format(r)
    if p_value >= 0.05:
        p_value_s = "(p > 0.05)"
    else:
        if p_value < 0.01:
            p_value_s = "(p < 0.01)"
        else:
            p_value_s = "(p < 0.05)"

    # number of valid values
    N_s = "N = " + str(n_colocations)

    # mean
    mean_x = np.mean(x_data)
    mean_y = np.mean(y_data)

    # standard deviation
    std_x  = np.std(x_data, ddof=1)
    std_y  = np.std(y_data, ddof=1)

    # mean and standard deviation
    x_mean_std_s = "x: " + str(round(mean_x, 3)) \
            + r" $\pm$ " + str(round(abs(std_x), 2))
    y_mean_std_s = "y: " + str(round(mean_y, 3)) \
            + r" $\pm$ " + str(round(abs(std_y), 2))

    # mean bias
    mb = ms.MB(x_data, y_data)
    mb_s = 'MB = {:.2f}'.format(mb)

    # normalized mean bias
    nmb = ms.NMB(x_data, y_data)
    nmb_s = 'NMB = {:.1f}%'.format(nmb*100)

    # mean squared error
    mse = ms.MSE(x_data, y_data)
    mse_s = 'MSE = {:.2f}'.format(mse)

    # mean absolute error
    mae = ms.MAE(x_data, y_data)
    mae_s = 'MAE = {:.3f}'.format(mae)

    # root mean squared error
    rmse = ms.RMSE(x_data, y_data)
    rmse_s = 'RMSE = {:.2f}'.format(rmse)

    # normalized centralized root mean squared error
    ncrmse = ms.NCRMSE(x_data, y_data)
    ncrmse_s = 'NCRMSE = {:.2f}'.format(ncrmse)

    # Normalized mean squared error. 
    # (normalized by mean(obs) * mean(est))
    nmse_mean = ms.NMSE_mean(x_data, y_data)
    nmse_mean_s = 'NMSE = {:.2f}'.format(nmse_mean)

    # Normalized mean squared error. 
    # (normalized by variance(obs))
    nmse_var = ms.NMSE_var(x_data, y_data)
    nmse_var_s = 'NMSE = {:.2f}'.format(nmse_var)

    # linear regression equation    
    sign = " + "
    if intercept < 0:
        sign = " - "
    linear_eq_s = "y = " + str(round(slope, 2)) + "x" + sign \
            + str(round(abs(intercept), 2))


    # statistic label dictornary
    stat_label = {
            'R': r_s + ' ' + p_value_s,
            'x_mean_std': x_mean_std_s,
            'y_mean_std': y_mean_std_s,
            'rmse': rmse_s,
            'ncrmse': ncrmse_s,
            'mse': mse_s,
            'mae': mae_s,
            'N': N_s,
            'mb': mb_s,
            'nmb': nmb_s,
            'nmse_mean': nmse_mean_s,
            'nmse_var': nmse_var_s,
            'linear_eq': linear_eq_s
            }
    label_list = [
            'N',
            'R',
            'linear_eq',
            'x_mean_std',
            'y_mean_std',
            'rmse',
            'ncrmse',
            'mse',
            'mae',
            'mb',
            'nmb',
            'nmse_mean',
            'nmse_var'
            ]
    print(' - scatter: statistics')
    for label in label_list:
        print(stat_label[label])

   
    # create scatter plot
    paths = ax.scatter(x_data, y_data, **kwargs)
    
    # add slope line
    min_x = np.nanmin(x_data)
    max_x = np.nanmax(x_data)
    min_y = np.nanmin(y_data)
    max_y = np.nanmax(y_data)
    x = np.array((min_x-3*np.absolute(min_x),max_x+3*np.absolute(max_x)))
    y = (slope * x) + intercept
    if regress_line:
        ax.plot(x, y, '-', color='black', linewidth=1.6)
    if one2one_line:
        ax.plot(x, x, '--', color='black', linewidth=0.8)
    
    # label_ul
    if label_ul is not None:
        label = ''
        for var_n in label_ul:
            var_v = stat_label[var_n]
            label = label + var_v + '\n'
        posXY0      = (0, 1)
        posXY_text0 = (5, -5)
        ax.annotate(label, xy=posXY0, xytext=posXY_text0, va='top', 
                xycoords='axes fraction', textcoords='offset points',
                color=text_color)

    # label_lr
    if label_lr is not None:
        label = ''
        for var_n in label_lr:
            var_v = stat_label[var_n]
            label = label + var_v + '\n'
        posXY0      = (1, 0)
        posXY_text0 = (-5, -5)
        ax.annotate(label, xy=posXY0, xytext=posXY_text0, 
                va='bottom',  ha='right',
                xycoords='axes fraction', textcoords='offset points',
                color=text_color)
    
    return paths, slope, intercept
#
#------------------------------------------------------------------------------
#
def plot_ee(ax, x, y, 
        high_slope=0.1, high_intercept=0.05,
        low_slope=0.1, low_intercept=0.05,
        vmin=0.0, vmax=5.0):
    """ Plot expected error envelop.
    +(high_intercept + high_slope * AOD)
    -(low_intercept  + low_slope  * AOD)
    (Yi Wang, 03/17/2021)
    """

    # expected error line
    xx = np.array([vmin, vmax])
    h_yy = xx + high_intercept + high_slope * xx
    l_yy = xx - low_intercept - low_slope * xx
    ax.plot(xx, h_yy, 'k--')
    ax.plot(xx, l_yy, 'k--')

    # stattistics according to expected error
    high_bound = x + high_intercept + high_slope * x
    high_result = (y-high_bound) > 0.0
    high_count = len(x[high_result])
    low_bound = x - low_intercept - low_slope * x
    low_result = (y-low_bound) < 0.0
    low_count = len(x[low_result])
    N = len(x)
    center_count = N - high_count - low_count
    if N > 0:
        high_ratio   = high_count   / float(N)
        low_ratio    = low_count    / float(N)
        center_ratio = center_count / float(N)
        label = 'within EE={:4.1f}%\n'.format(center_ratio*100) + \
                'above EE={:4.1f}%\n'.format(high_ratio*100) + \
                'below EE={:4.1f}%'.format(low_ratio*100)
        posXY0      = (1, 0)
        posXY_text0 = (-5, 5)
        ax.annotate(label, xy=posXY0, xytext=posXY_text0,
                va='bottom',  ha='right',
                xycoords='axes fraction', textcoords='offset points',
                color='k')
#
#------------------------------------------------------------------------------
#
