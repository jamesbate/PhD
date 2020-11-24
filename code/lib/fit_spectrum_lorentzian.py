"""This function fits lorentzians to spectral peaks. 
Maybe I should generalise multiple peak fitting?
"""
##-------------------------------PREAMBLE-----------------------------------##
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq, minimize
import numpy as np
#external packages
##-------------------------------PARAMETERS-----------------------------------##
def lorentzian(x, x_0, a, gam, cutoff = None):
    @np.vectorize
    def square(x, c):
        if a * gam**2 / ( gam**2 + ( x - x_0 )**2) > c:
            return 0
        return 1
    if cutoff is not None:
        return  square(x, cutoff) * a * gam**2 / ( gam**2 + ( x - x_0 )**2)
    return  a * gam**2 / ( gam**2 + ( x - x_0 )**2)

##-------------------------------FIT-----------------------------------##

def multi_lorentzian(x, params, cutoff = None):

    offset = params[0]
    rest_params = params[1:]
    return offset + sum([lorentzian(x, *rest_params[i :i+3], cutoff = cutoff) for i in range(0, len(rest_params),3)])

def res_multi_lorentzian(params, x_data, y_data):
    diff = [ multi_lorentzian( x, params ) - y for x, y in zip( x_data, y_data )]
    return diff

def fit_spectrum_lorentzian(freq, probs, dprobs,initial_parameters, labels = None):
    """
    Capable of fitting/plotting for multiple spectra, so inputs
    should be array like
    """
    plt.rcParams.update({'font.size': 16})
    colours = ['b','m','c','r','tab:orange', 'tab:pink']
    fig, ax = plt.subplots(1	,1, constrained_layout=True, figsize=(18, 9))
    plt.grid()
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Excitation Probability')
    plt.title('729 Spectroscopy')
    #figure settings

    nplots = freq.shape[0]
    #number of separate plots

    fit_values = []
    fitting_params = []
    loop = 0

    for f,p,d,i_p in zip(freq,probs,dprobs, initial_parameters):

        ax.scatter(f, p, color = 'k', s = 5)
        #plot data

        ax.errorbar(f, p, d, ls = 'none', c = 'k', capsize = 3)
        #errorbars

        #find guesses for fitting params
        nmax = np.argmax(probs)
        nmin = np.argmin(probs)
        maxfreq = f[nmax]
        minfreq = f[nmin]

        init_params = [min(p)]
        #location of maximum
        f_0 = f[nmax]

        init_params += i_p


        #fit lorenztians
        popt, cov_x,_,_,_ = leastsq( res_multi_lorentzian, init_params, args=(f, p), full_output = 1)
        #record params
        fitting_params.append([popt, cov_x])

        #plot fit
        fit_domain = np.linspace(f[0], f[-1], 1000)
        if labels is not None:
            ax.plot(fit_domain, multi_lorentzian(fit_domain, popt, cutoff = 30), color = colours[loop], label = labels[loop], linewidth = 4)
        else:
            ax.plot(fit_domain, multi_lorentzian(fit_domain, popt, cutoff = 30), color = colours[loop], linewidth = 4)
        loop += 1
    if labels is not None:
        plt.legend()

    return fitting_params
