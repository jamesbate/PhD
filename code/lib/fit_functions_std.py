"""Since parameters are passed by a dictionary, 
default values must be implemented manually 
"""
##-------------------------------PREAMBLE-----------------------------------##

import numpy as np 
import matplotlib.pyplot as plt  
##---------------------------FITTING FUNCTIONS-------------------------------##

def fit_sinusoid(x, parameters_dict):
    """function: A*np.cos(B*x + C) 
    default A: 1 
    default C: 0
    """
    #necessary
    B = parameters_dict['B']

    #defaulted
    if 'A' in parameters_dict:
        A = parameters_dict['A']
    else: 
        A = 1
    if 'C' in parameters_dict:
        C = parameters_dict['C']
    else: 
        C = 0

    #function
    return A*np.cos(B*x + C) 

def fit_poissonian(x,mean,A):
    """function: A*(mean**x*np.exp(-mean))/(factorial(x))
    default A: 1
    """
    #necessary
    mean = parameters_dict['mean']

    #defaulted
    if 'A' in parameters_dict:
        A = parameters_dict['A']
    else: 
        A = 1
    #function
	return A*(mean**x*np.exp(-mean))/(factorial(x))
