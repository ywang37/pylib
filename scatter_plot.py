import matplotlib.pyplot as plt
import numpy as np

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

def scatter(ax, in_x_data, in_y_data, label_p = 'upper left', \
        regress_line=True, one2one_line=True, **kwargs):
    """ Function to create scatter plots.
    """
    
    # copy and flatten data
    x_data = in_x_data.flatten()
    y_data = in_y_data.flatten()

    # get a mask with those elements posible to compare (non-nan values)
    mask = np.logical_and(np.logical_not(np.isnan(x_data)), 
                          np.logical_not(np.isnan(y_data)))
    n_colocations = len(mask[mask==True])
    x_data = x_data[mask]
    y_data = y_data[mask]
    
    # liner regression
    slope, intercept, correlation, p_value_slope, std_error = \
            stats.linregress(x_data, y_data)

    # Calculates a Pearson correlation coefficient 
    # and the p-value for testing non-correlation
    r, p_value = stats.pearsonr(x_data, y_data)

    rmse   = np.std(y_data-x_data, ddof=0)
    mean_x = np.mean(x_data)
    mean_y = np.mean(y_data)
    std_x  = np.std(x_data, ddof=1)
    std_y  = np.std(y_data, ddof=1)
    
    # create scatter plot
    paths = ax.scatter(x_data, y_data, **kwargs)
    
    # add slope line
    min_x = np.nanmin(x_data)
    max_x = np.nanmax(x_data)
    min_y = np.nanmin(y_data)
    max_y = np.nanmax(y_data)
    x = np.array((min_x-np.absolute(min_x),max_x+np.absolute(max_x)))
    y = (slope * x) + intercept
    if regress_line:
        ax.plot(x, y, '-', color='black', linewidth=1.6)
    if one2one_line:
        ax.plot(x, x, '--', color='black', linewidth=0.8)
    
    # create strings for equations in the plot
    correlation_string = "R = {:.2f}".format(r)
    
    sign = " + "
    if intercept < 0:
        sign = " - "
        
    lineal_eq = "y = " + str(round(slope, 2)) + "x" + sign \
            + str(round(abs(intercept), 2))
    rmse_coef = "RMSE = " + str(round(rmse,2))

    if p_value >= 0.05:
        p_value_s = "(p > 0.05)"
    else:
        if p_value < 0.01:
            p_value_s = "(p < 0.01)"
        else:
            p_value_s = "(p < 0.05)"

    n_collocations = "N = " + str(n_colocations)
    x_mean_std = "x: " + str(round(mean_x, 3)) \
            + " $\pm$ " + str(round(abs(std_x), 2))
    y_mean_std = "y: " + str(round(mean_y, 3)) \
            + " $\pm$ " + str(round(abs(std_y), 2))

    print(correlation_string + ' ' + p_value_s)

    if (label_p == 'upper left'):
        equations0 = \
                correlation_string + ' ' + p_value_s + '\n' + \
                x_mean_std + '\n' + \
                y_mean_std + '\n' + \
                rmse_coef  + '\n' + \
                n_collocations

        # upper left
        posXY0      = (0, 1)
        posXY_text0 = (5, -5)
        ax.annotate(equations0, xy=posXY0, xytext=posXY_text0, va='top', \
                xycoords='axes fraction', textcoords='offset points')
    elif (label_p == 'lower right'):
        equations0 = n_collocations + '\n' + \
                rmse_coef  + '\n' + \
                x_mean_std + '\n' + \
                y_mean_std # + '\n' + \
                #lineal_eq  + '\n' + \
                #correlation_string + ' ' + p_value_s

        # lower right
        posXY0      = (1, 0)
        posXY_text0 = (-5, 5)
        ax.annotate(equations0, xy=posXY0, xytext=posXY_text0, 
                va='bottom', ha='right',
                xycoords='axes fraction', textcoords='offset points')

    elif (label_p == 'positive'):
        equations1 = correlation_string + ' ' + p_value_s + '\n'+ \
                lineal_eq  + '\n' + \
                rmse_coef  + '\n'
        equations2 = n_collocations + '\n' + \
                x_mean_std + '\n' + \
                y_mean_std
            
        posXY1      = (0, 1)
        posXY_text1 = (5, -5)
        ax.annotate(equations1, xy=posXY1, xytext=posXY_text1, va='top', \
                xycoords='axes fraction', textcoords='offset points')
    
        posXY2      = (1, 0)
        posXY_text2 = (-5, 5)
        ax.annotate(equations2, xy=posXY2, xytext=posXY_text2, 
                va='bottom', ha='right', 
                xycoords='axes fraction', textcoords='offset points')
    else:
        '!!! scatter: label_p error !!!'
        exit()
    
    return paths, slope, intercept
