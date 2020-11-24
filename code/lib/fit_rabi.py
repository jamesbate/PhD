"""This function plots rabi frequencies with error bars, fitting function
uses full sideband frequency, not LD approximation
"""
##-------------------------------PREAMBLE-----------------------------------##
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq, minimize
import numpy as np
from scipy.special import eval_genlaguerre, factorial
#external packages

import math
#stdlib imports

from .physical_constants import PhysicalConstants
from .statistical_physics import therm_prob
#internal imports
##-------------------------------PARAMETERS-----------------------------------##

#Let our ion be in some thermal state

def full_rabi_dicke(n, m, dicke_factor, Rabi):
	#returns full coupling strength for any sideband, not LD regime

	#cannot couple below motional 0 state, introduces asymmetry
	if n + m < 0:
		return 0
		#RSB cannot be driven below ground thermal state
	if m > 0:
		n_2 = n + m
		n_1 = n
	else:
		n_2 = n
		n_1 = n + m
	#convention here is that n_2 is bigger and n_1 is smaller
	pf = 1/partial_factorial(n_2, n_1)
	return abs(Rabi*np.exp(-dicke_factor**2/2)*dicke_factor**(abs(m))*eval_genlaguerre(n_1, abs(m), dicke_factor**2)*(pf)**(1/2))

detuned_rabi_factor = lambda rabi_dicke, detuning: rabi_dicke/np.sqrt(rabi_dicke**2 + detuning**2)
pc = PhysicalConstants()

def full_Rabi_wrad(t, p, sum_limit = 100, p_fix = [None,None,None,None,None, None, None], m = [0,0]):
	#m represents sideband, DO NOT USE ABOVE 3. Assumes thermal distribution
	# if m > 3:
	# 	raise ValueError('fix full_rabi_dicke before m > 3')

    #unpack parameters 
    dicke_factor_r = parameters_dict['dicke_factor_r']
    dicke_factor_a = parameters_dict['dicke_factor_a']
    n_av_a = parameters_dict['n_av_a']
    n_av_r = parameters_dict['n_av_r']
    detuning = parameters_dict['detuning']
    Rabi = parameters_dict['Rabi']
    phase = parameters_dict['phase']
    amplitude = parameters_dict['amplitude']



	#and now calculate function
    ret = 0
    for n_1 in range(0,sum_limit):
    	for n_2 in range(0,sum_limit):
            dicke = full_rabi_dicke(n_1, m[0], dicke_factor_a, Rabi)
            dicke *= full_rabi_dicke(n_2, m[1], dicke_factor_r, Rabi)**2/Rabi**2

            if dicke == 0:
                continue
            det_fac = detuned_rabi_factor(dicke, detuning)
            ret += amplitude*therm_prob(n_1, n_av_a)*therm_prob(n_2, n_av_r)*det_fac**2*np.sin(phase + 0.5*t*dicke/det_fac)**2
    return np.array([ret]).T

def full_Rabi(t, parameters_dict, **kwargs):
    """Say stuff about this function 
    Note: m represents sideband, DO NOT USE ABOVE 3. Assumes thermal distribution
    """
    m = 0 
    sum_limit = 100
    #defaulting parameters: 
    if 'm' in kwargs.keys():
        m = kwargs['m']
    if 'sum_limit' in kwargs.keys():
        sum_limit = kwargs['sum_limit']

    if m > 3:
        raise ValueError('fix full_rabi_dicke before m > 3')

    #unpack parameters 
    dicke_factor = parameters_dict['dicke_factor']
    n_av = parameters_dict['n_av']
    detuning = parameters_dict['detuning']
    Rabi = parameters_dict['Rabi']
    phase = parameters_dict['phase']
    amplitude = parameters_dict['amplitude']

    ret = 0
    for n in range(0,sum_limit):
        if full_rabi_dicke(n, m, dicke_factor, Rabi) == 0:
            continue
        det_fac = detuned_rabi_factor(full_rabi_dicke(n, m, dicke_factor, Rabi), detuning)
        ret += amplitude*therm_prob(n, n_av)*det_fac**2*np.sin(phase + 0.5*t*full_rabi_dicke(n, m, dicke_factor, Rabi)/det_fac)**2

    return np.array([ret]).T



def partial_factorial(n,m):
    #not including m
    tot = 1
    for i in range(0, n-m):
        tot *= n-i
    return tot

def COM_freq(lamb_dicke):
	#assuming our experimental configuration
	#factor of 2 in denom from angle
	k = ((2*np.pi)/729e-9)
	return (pc.h_bar*k**2/(2*np.pi*2*2*pc.m_calc*lamb_dicke**2))

