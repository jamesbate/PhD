"""Defines functions of excitation of carrier and first order sidebands in the
 Lamb-~Dicke regime, taking up to first order in n.
"""
##-------------------------------PREAMBLE-----------------------------------##
import numpy as np
import matplotlib.pyplot as plt
#external dependencies

##-------------------------------PARAMETERS-----------------------------------##

#when to cap the sum over motional modes
sum_limit = 500

##--------------------------FUNCTION DEFINITIONS------------------------------##

def rabi_flops_LD_red(t, rabi, eta, n_av):
    r = 0
    for n in range(sum_limit):
        r += ((n_av)/(n_av + 1))**(n+1)*(np.sin(rabi*eta*t*np.sqrt(n + 1)/2)**2)/(n_av + 1)
    return r

def rabi_flops_LD_blue(t, rabi, eta, n_av):
    r = 0
    for n in range(sum_limit):
        r += ((n_av)/(n_av + 1))**(n)*(np.sin(rabi*eta*t*np.sqrt(n + 1)/2)**2)/(n_av + 1)
    return r

def rabi_flops_LD_carrier(t, rabi, eta, n_av):
    r = 0
    for n in range(sum_limit):
        r += ((n_av)/(n_av + 1))**(n)*(np.sin(rabi*t*(1 - n*eta**2)/2)**2)/(n_av + 1)
    return r
