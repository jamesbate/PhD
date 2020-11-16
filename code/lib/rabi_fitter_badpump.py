"""This function plots rabi frequencies with error bars and with a fitted
damped sinusoid.
"""
##-------------------------------PREAMBLE-----------------------------------##
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq, minimize
import numpy as np
from .physical_constants import PhysicalConstants
from scipy.special import eval_genlaguerre, factorial
import math 
#external packages
##-------------------------------PARAMETERS-----------------------------------##

therm_prob = lambda n, n_av: (n_av/(n_av + 1))**n/(n_av + 1)
#Let our ion be in some thermal state
Rabi_dicke = lambda n, dicke_factor, Rabi: Rabi*(1 - dicke_factor**2*n)
#Second order carrier coupling strength in LD regime, NOT USED NOW

def full_rabi_dicke(n, m, dicke_factor, Rabi):
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
#full coupling strength for any sideband, not LD regime


detuned_rabi_factor = lambda rabi_dicke, detuning: rabi_dicke/np.sqrt(rabi_dicke**2 + detuning**2)
#might not use this

pc = PhysicalConstants()

def full_Rabi(t, p, sum_limit = 100, p_fix = [None,None,None,None,None, None], m = 0):
    #m represents sideband, DO NOT USE ABOVE 3
    if m > 3: 
        raise ValueError('fix full_rabi_dicke before m > 3')
    p_temp = [None,None,None,None,None, None]
    p_count = 0
    for n,i in enumerate(p_fix):
        if i is not None:
            p_temp[n] = i 
        else:
            p_temp[n] = p[p_count]
            
            p_count += 1

    [dicke_factor, n_av, detuning, Rabi, phase, A] = p_temp

    if n_av < 0:
        return 0
        #unphysical

    ret = 0

    for n in range(0,sum_limit):
        if full_rabi_dicke(n, m, dicke_factor, Rabi) == 0:
            continue
        det_fac = detuned_rabi_factor(full_rabi_dicke(n, m, dicke_factor, Rabi), detuning)
        ret += A*therm_prob(n, n_av)*det_fac**2*np.sin(phase + 0.5*t*full_rabi_dicke(n, m, dicke_factor, Rabi)/det_fac)**2
    return ret

#physical constraints for minimisation
##-------------------------FUNCTION DEFINITION--------------------------------##

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


def rabi_fitter_badpump(phase, probs, errors, initial_fit, labels = None, p_fix = [None,None,None,None,None, None],m = 0, ax = None, color = None):
	"""This function takes raw rabi data, calculates error bars, fits a dambed
	sinusoid, and creates a plot figure. It returns fitting parameters.
	Parameters
	----------
	phase:
		array-like, become x axis of plot. Typically will be the phase of the
		qubit laser i.e. how long the pulse lasted for

	probs:
		2d numpy array, each column is taken to be a dataset. (for one set of
		data, must in form [[data1],[data2],[data3],...]).
	labels:
		array-like, at least as long as the number of different rabi datasets,
		will be indexed for each dataset and passed as colour argument for plot.
	initial_fit:
		2d array_like, initial try for scipy curve_fit for each dataset, suggested
	Returns
	-------
	fit_values:
		list, parameters for each fit
	"""
	#plt.rcParams.update({'font.size': 16})
	colours = ['b','m','c','r','tab:orange', 'tab:pink']
	if ax is None:
		fig, ax = plt.subplots(1	,1, constrained_layout=True, figsize=(18, 9))

		plt.grid()
		#figure settings

	fit_values = []
	(c,n) = probs.shape

	if labels is None:
		labels = ['p'+str(i) for i in range(n)]
	#default legend labels if none provided (not really recommended)

	for i in range(n):
		#loop over datasets

		def func_to_minimise(p):
			return np.sum((full_Rabi(phase, p, p_fix = p_fix[i], m = m) - probs[:,i])**2)
		#define least squares metric to minimize
		if color is None:
			ax.scatter(phase, probs[:,i], c = colours[i], label = labels[i])

			#calculate projection noise
			ax.errorbar(phase, probs[:,i], errors[:,i], ls = 'none', c = colours[i], capsize = 3)
			#plot errorbars
		else: 
			ax.scatter(phase, probs[:,i], c = color[i], label = labels[i])

			#calculate projection noise
			ax.errorbar(phase, probs[:,i], errors[:,i], ls = 'none', c = color[i], capsize = 3)
			#plot errorbars
		opt_params = [initial_fit[i][j] for j in range(len(initial_fit[i])) if p_fix[i][j] is None]

		res = minimize(func_to_minimise, x0 = opt_params)
		#fit function with errors

		ftol = 2.220446049250313e-09
		#default tolerance 
		tmp_i = np.zeros(len(res.x))
		uncertainty = np.zeros(len(res.x))
		for k in range(len(res.x)):
			tmp_i[k] = 1.0
			hess_inv_i = res.hess_inv.dot(tmp_i)[k]
			uncertainty[k] = np.sqrt(max(1, abs(res.fun)) * ftol * hess_inv_i)
			tmp_i[k] = 0.0

		#now append full set of parameters again 
		final_result = [] 
		res_count = 0
		for j in range(len(initial_fit[i])):
			if p_fix[i][j] is None:
				final_result.append([res.x[res_count], uncertainty, full_Rabi(phase, res.x, p_fix = p_fix[i], m = m)])
				res_count += 1 
			else: 
				final_result.append([p_fix[i][j], 0, full_Rabi(phase, res.x, p_fix = p_fix[i], m = m)])
		fit_label = "Rabi: " + str(round(final_result[3][0]/2,2)) + "MHz\nPhonons: " + str(round(final_result[1][0],2)) + "\nAmplitude: " + str(round(final_result[5][0],2))

		fitdomain = np.linspace(phase[0], phase[-1], 1000)
		if color is None:
			ax.plot(fitdomain, full_Rabi(fitdomain, res.x, p_fix = p_fix[i], m = m), color = colours[i], label = fit_label)
		else:
			ax.plot(fitdomain, full_Rabi(fitdomain, res.x, p_fix = p_fix[i], m = m), color = color[i], label = fit_label)
	
		#ax.plot(fitdomain, full_Rabi(fitdomain, initial_fit[i]), c = 'k')
		fit_values.append(final_result)
		print("COM Frequency (MHz): ",COM_freq(final_result[0][0])/1e6)

	ax.legend()
	ax.set_title('Rabi Flops')
	ax.set_xlabel('Pulse Length (us)')
	ax.set_ylabel('Mean Excitation')

	return fit_values
