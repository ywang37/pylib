import numpy as np
from scipy.optimize import curve_fit

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
        return a + b * x + c11 * np.sin(2.0*np.pi*1*x) + c12 * np.cos(2.0*np.pi*1*x)

    @staticmethod
    def model_2(x, a, b, c11, c12, c21, c22):
        """
        """
        return a + b * x + c11 * np.sin(2.0*np.pi*1*x) + c12 * np.cos(2.0*np.pi*1*x) + \
                           c21 * np.sin(2.0*np.pi*2*x) + c22 * np.cos(2.0*np.pi*2*x)
                           
    @staticmethod
    def model_3(x, a, b, c11, c12, c21, c22, c31, c32):
        """
        """
        return a + b * x + c11 * np.sin(2.0*np.pi*1*x) + c12 * np.cos(2.0*np.pi*1*x) + \
                           c21 * np.sin(2.0*np.pi*2*x) + c22 * np.cos(2.0*np.pi*2*x) + \
                           c31 * np.sin(2.0*np.pi*3*x) + c32 * np.cos(2.0*np.pi*3*x)

    @staticmethod
    def model_4(x, a, b, c11, c12, c21, c22, c31, c32, c41, c42):
        """
        """
        return a + b * x + c11 * np.sin(2.0*np.pi*1*x) + c12 * np.cos(2.0*np.pi*1*x) + \
                           c21 * np.sin(2.0*np.pi*2*x) + c22 * np.cos(2.0*np.pi*2*x) + \
                           c31 * np.sin(2.0*np.pi*3*x) + c32 * np.cos(2.0*np.pi*3*x) + \
                           c41 * np.sin(2.0*np.pi*4*x) + c42 * np.cos(2.0*np.pi*4*x)

    def analysis(self, x, y):
        """
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
        self.r1        = np.corrcoef(self.noise[0:n_noise-1], self.noise[1:n_noise])[0,1]
        self.noise_std = np.std(self.noise)

        # 
        self.epsilon_std = np.sqrt(self.noise_std**2 * (1-self.r1**2))

        # standard deviation of trend
        self.trend_std = self.noise_std / ((len(x)/12.0)**1.5) * np.sqrt((1.0+self.r1)/(1.0-self.r1))


        
















        

