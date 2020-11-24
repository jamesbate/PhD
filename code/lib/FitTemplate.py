"""After seeing how great the lmfit package, I was inspired to create my own
object using it. This acts as a fitting template. 
"""
##-------------------------------PREAMBLE-----------------------------------##
import numpy as np 
import matplotlib.pyplot as plt 
from lmfit import minimize, Parameters, fit_report 
import logging 

##-------------------------------CLASS DEFINITION-----------------------------------##

class FitTemplate(): 
    def __init__(self, fit_function, log_dir = None):
        self.fit_function = fit_function 
        self.parameters = Parameters()
        self.fit_result = None

        #setup logging. warning level is standard and is sent to stdout. info is requested by log_dir argument,
        #and is printed to log file
        if log_dir is not None: 
            logging.basicConfig(filename=log_dir +'log.log', level=logging.INFO)
        else:
            logging.basicConfig(level=logging.CRITICAL)
    

    def residuals_wrapper(self, parameters, x, data,weights,**kwargs):
        model_values = self.fit_function(x, parameters.valuesdict(), **kwargs)
        return ((model_values - data)*weights)**2
        
    def do_minimisation(self, x, data, weights = 1, **kwargs):
        self.fit_result = minimize(self.residuals_wrapper, self.parameters, args = (x, data, weights), kws = kwargs)
        logging.info('Fit Result')
        logging.info('==========')
        return self.fit_result

    def get_opt_parameters(self):
        if self.fit_result is None: 
            raise ValueError("No fit result! Do a fit before asking for")
        return self.fit_result.params.valuesdict()

    def print_parameters(self):
        self.parameters.pretty_print() 
    
    def print_fit_result(self):
        logging.info((fit_report(self.fit_result)))
        print(fit_report(self.fit_result))

    def plot_fit(self, x, y, xlabel = None, ylabel = None, title = None, errorbars = None, label = None, ax = None, c = None, colour_index = None, **kwargs): 

        if ax is None:
            _, ax = plt.subplots(1	,1, constrained_layout=True, figsize=(18, 9))
        plt.rcParams.update({'font.size': 16})       
        colours = ['b','m','c','r','tab:orange', 'tab:pink']

        #decide colour 
        if c is not None: 
            color = c 
        elif colour_index is not None: 
            color = colours[colour_index]
        else: 
            color = colours[0]

        #scatter plot
        ax.scatter(x, y, color = color)
        #plot errors
        if errorbars is not None:
            ax.errorbar(x, y, errorbars, ls = 'none', c = color, capsize = 3)
        #plot model
        fitdomain = np.linspace(x[0], x[-1], 1000)	
        ax.plot(fitdomain, self.fit_function(fitdomain, self.fit_result.params.valuesdict(), **kwargs), c = color, label = label)
        plt.legend()
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        plt.grid()
        return ax 
        
        		
        
        
        





