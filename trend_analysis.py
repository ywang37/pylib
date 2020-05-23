from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from mylib.cartopy_plot import cartopy_plot
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.io import read_nc

#
#------------------------------------------------------------------------------
#
class trend_analysis():
    """
    Trend analysis
    """

    def __init__(self, fit_model_index=0):
        self.fit_model_index = fit_model_index
        self.fit_model       = None
        self.popt            = None
        self.pcov            = None
        self.y_fit           = None
        self.y_trend         = None
        self.y_period        = None
        self.noise           = None
        self.r1              = None
        self.noise_std       = None
        self.epsilon_std     = None
        self.trend_std       = None

    @staticmethod
    def model_0(x, a, b):
        """
        """
        return a + b * x

    @staticmethod
    def model_1(x, a, b, c11, c12):
        """
        """
        return a + b * x + c11 * np.sin(2.0*np.pi*1*x) + \
                           c12 * np.cos(2.0*np.pi*1*x)

    @staticmethod
    def model_2(x, a, b, c11, c12, c21, c22):
        """
        """
        return a + b * x + c11 * np.sin(2.0*np.pi*1*x) + \
                           c12 * np.cos(2.0*np.pi*1*x) + \
                           c21 * np.sin(2.0*np.pi*2*x) + \
                           c22 * np.cos(2.0*np.pi*2*x)
                           
    @staticmethod
    def model_3(x, a, b, c11, c12, c21, c22, c31, c32):
        """
        """
        return a + b * x + c11 * np.sin(2.0*np.pi*1*x) + \
                           c12 * np.cos(2.0*np.pi*1*x) + \
                           c21 * np.sin(2.0*np.pi*2*x) + \
                           c22 * np.cos(2.0*np.pi*2*x) + \
                           c31 * np.sin(2.0*np.pi*3*x) + \
                           c32 * np.cos(2.0*np.pi*3*x)

    @staticmethod
    def model_4(x, a, b, c11, c12, c21, c22, c31, c32, c41, c42):
        """
        """
        return a + b * x + c11 * np.sin(2.0*np.pi*1*x) + \
                           c12 * np.cos(2.0*np.pi*1*x) + \
                           c21 * np.sin(2.0*np.pi*2*x) + \
                           c22 * np.cos(2.0*np.pi*2*x) + \
                           c31 * np.sin(2.0*np.pi*3*x) + \
                           c32 * np.cos(2.0*np.pi*3*x) + \
                           c41 * np.sin(2.0*np.pi*4*x) + \
                           c42 * np.cos(2.0*np.pi*4*x)

    def analysis(self, x, y):
        """ for monthly data.
        """
        
        if self.fit_model_index == 0:
            self.fit_model = self.model_0
        elif self.fit_model_index == 1:
            self.fit_model = self.model_1
        elif self.fit_model_index == 2:
            self.fit_model = self.model_2
        elif self.fit_model_index == 3:
            self.fit_model = self.model_3
        elif self.fit_model_index == 4:
            self.fit_model = self.model_4
        else:
            print('ERROR')
            exit()
        
        # get model parameters
        popt, pcov = curve_fit(self.fit_model, x, y)
        self.popt = popt
        self.pcov = pcov

        # reconstruct model
        if self.fit_model_index == 0:
            self.y_fit = self.fit_model(x, self.popt[0], self.popt[1])
        if self.fit_model_index == 1:
            self.y_fit = self.fit_model(x, self.popt[0], self.popt[1],\
                                           self.popt[2], self.popt[3])
        if self.fit_model_index == 2:
            self.y_fit = self.fit_model(x, self.popt[0], self.popt[1],\
                                           self.popt[2], self.popt[3],\
                                           self.popt[4], self.popt[5])
        if self.fit_model_index == 3:
            self.y_fit = self.fit_model(x, self.popt[0], self.popt[1],\
                                           self.popt[2], self.popt[3],\
                                           self.popt[4], self.popt[5],\
                                           self.popt[6], self.popt[7])
        if self.fit_model_index == 4:
            self.y_fit = self.fit_model(x, self.popt[0], self.popt[1],\
                                           self.popt[2], self.popt[3],\
                                           self.popt[4], self.popt[5],\
                                           self.popt[6], self.popt[7],\
                                           self.popt[8], self.popt[9])

        # linear trend
        self.y_trend  = self.model_0(x, self.popt[0], self.popt[1])
        self.y_period = self.y_fit - self.y_trend

        # calculate noise
        self.noise = y - self.y_fit

        # noise autocorrelation coefficient and standard deviation
        n_noise        = len(self.noise)
        self.r1        = np.corrcoef(self.noise[0:n_noise-1], \
                                     self.noise[1:n_noise])[0,1]
        self.noise_std = np.std(self.noise)

        # 
        self.epsilon_std = np.sqrt(self.noise_std**2 * (1-self.r1**2))

        # standard deviation of trend
        self.trend_std = self.noise_std / \
                ((len(x)/12.0)**1.5) * np.sqrt((1.0+self.r1)/(1.0-self.r1))

    def analysis_yearly(self, y):
        """ for year data
        (ywang, 05/22/20)
        """

        # no seasonal variablity
        self.fit_model = self.model_0

        # get model parameters
        x = np.array(range(len(y)), dtype=float)
        popt, pcov = curve_fit(self.fit_model, x, y)
        self.popt = popt
        self.pcov = pcov

        # reconstruct model
        self.y_fit = self.fit_model(x, self.popt[0], self.popt[1])

        # linear trend
        self.y_trend  = self.model_0(x, self.popt[0], self.popt[1])
        self.y_period = self.y_fit - self.y_trend

        # calculate noise
        self.noise = y - self.y_fit

        # noise autocorrelation coefficient and standard deviation
        n_noise        = len(self.noise)
        self.r1        = np.corrcoef(self.noise[0:n_noise-1], \
                                     self.noise[1:n_noise])[0,1]
        self.noise_std = np.std(self.noise)

        # 
        self.epsilon_std = np.sqrt(self.noise_std**2 * (1-self.r1**2))

        # standard deviation of trend
        self.trend_std = self.noise_std / \
                (len(x)**1.5) * np.sqrt((1.0+self.r1)/(1.0-self.r1))
#
#------------------------------------------------------------------------------
#
def plot_trend_map(filename, fig_dir, mean_flag=True, trend_flag=True, 
        sigma_flag=True, name='', 
        mean_vmin=None, mean_vmax=None,
        mean_units='',
        trend_vmin=None, trend_vmax=None, 
        trend_cmap=plt.get_cmap('seismic'), trend_units='',
        sigma_units='',
        ):
    """ Plot results from trend_analysis.
    (ywang, 05/21/20)

    Parameters
    ----------
    filename : str
        netCDF file of trend_analysis results.
    fig_dir : str
        Directory to save figures.
    mean_flag : bool
        Plot mean or not
    trend_flag : bool
        Plot trend or not
    sigma_flag : bool
        Plot trend standard deviation or not
    name : str
        Prefix of figure names.

    """

    # directory
    if fig_dir[-1] != '/':
        fig_dir = fig_dir + '/'

    # variables to be read
    varnames = ['Latitude_e', 'Longitude_e']
    if mean_flag:
        varnames.append('mean')
    if trend_flag:
        varnames.append('trend')
    if sigma_flag:
        varnames.append('trend_sigma')

    # read data
    data = read_nc(filename, varnames, verbose=True)

    # get latitude and longitude edges
    lat_e = data['Latitude_e']
    lon_e = data['Longitude_e']

    cbar_prop = {}
    cbar_prop['orientation'] = 'horizontal'

    # plot mean
    if mean_flag:

        mean = data['mean']

        mean_p = cartopy_plot(lon_e, lat_e, mean, cbar_prop=cbar_prop,
                vmin=mean_vmin, vmax=mean_vmax)
        mean_p['cb'].set_label(mean_units)

        fig_mean = fig_dir + name + '_mean.png'
        plt.savefig(fig_mean, format='png', dpi=300)

    # plot trend
    if trend_flag:

        trend = data['trend']

        trend_p = cartopy_plot(lon_e, lat_e, trend, cbar_prop=cbar_prop,
                vmin=trend_vmin, vmax=trend_vmax,
                cmap=trend_cmap)
        trend_p['cb'].set_label(trend_units)

        fig_trend = fig_dir + name + '_trend.png'
        plt.savefig(fig_trend, format='png', dpi=300)

    # plot trend_sigma
    if sigma_flag:

        sigma = data['trend_sigma']

        if not trend_flag:
            trend = data['trend']

        # plot sigma
        sigma_p = cartopy_plot(lon_e, lat_e, sigma, cbar_prop=cbar_prop)
        sigma_p['cb'].set_label(sigma_units)

        # save sigma plot
        fig_sigma = fig_dir + name + '_sigma.png'
        plt.savefig(fig_sigma, format='png', dpi=300)

        # trends that are significant
        trend_signi_flag = (np.absolute(trend) / sigma) >= 2.0
        trend_signi = \
                np.ma.masked_array(trend, np.logical_not(trend_signi_flag))

        # plot trends that are significant
        trend_signi_p = cartopy_plot(lon_e, lat_e, trend_signi, 
                cbar_prop=cbar_prop,
                vmin=trend_vmin, vmax=trend_vmax,
                cmap=trend_cmap)
        trend_signi_p['cb'].set_label(trend_units)


        # save trends that are significant plot
        fig_trend_signi = fig_dir + name + '_trend_signi.png'
        plt.savefig(fig_trend_signi, format='png', dpi=300)




        








#
#------------------------------------------------------------------------------
#
