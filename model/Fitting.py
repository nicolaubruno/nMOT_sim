# Libraries
import numpy as np

from scipy.optimize import curve_fit
from scipy.special import erf

#
class Fitting:
    #--
    # Attributes
    #--

    @property
    def params(self):
        return self._params

    @property
    def data(self):
        return self._data

    @property
    def func(self):
        return self._func
    
    #--
    # Methods
    #--

    #
    def __init__(self, X: list, Y: list, fit_func: str):
        # Set data
        self._data = {"X": X, "Y": Y}

        # Polynomial fitting
        #--
        poly_fit = False
        for i in range(1, 11):
            if fit_func == ("%ddeg_poly" % i):
                degree = i
                poly_fit = True
                break

        if poly_fit:
            self.poly(degree)
        #--

        # Invalid fitting function
        else: raise ValueError('Fitting function is not available')

    # Polynomial fitting
    def poly(self, degree):
        # Fit polynomial
        self._params = np.polyfit(self.data["X"], self.data["Y"], degree, full = False, cov = False)

        # Fitting function
        self._func = lambda x: np.vdot(self.params, np.array([x**n for n in range(degree, -1, -1)]))

    # Error function fitting
    def erf(self):
        # Error function
        f = lambda x, amp, mean, std: amp*(1 - (erf((t - mean) / np.sqrt(2 * std**2)) - erf((- mean) / np.sqrt(2 * std**2))) / 2)

        # Fitting
        max_X = np.max(self.data["X"])
        min_X = np.min(self.data["X"])
        self._params, covs = curve_fit(f, self.data["X"], self.data["Y"], bounds=([0, min_X, 0], [1, max_X, (max_X - min_X)]))

        self._func = lambda x: f(x, self.params[0], self.params[1], self.params[2])