import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from scipy.special import factorial
from .create_hist_plot import create_hist_plot

def PMT_fitter(data, labels = None, thresholds = None, maxim = 100, minim = 0, binsize = 1, ax = None):

	plt.rcParams.update({'font.size': 16})

	colours = ['b','m','c','r','tab:orange', 'tab:pink']

	fit_values = []

	binedges = np.arange(minim - binsize/2, maxim + binsize/2, binsize)
	domain = np.arange(minim,maxim, binsize)

	def poissonian(x,mean,A):
		return A*(mean**x*np.exp(-mean))/(factorial(x))

	if thresholds == None:
		nthresh = 1
		thresh_edges = [minim, maxim]
		guesses = [0,0.1]
	else:
		nthresh = len(thresholds) + 1
		thresh_edges = [minim] + thresholds + [maxim]
		guesses = [thresh_edges[i+1] - thresh_edges[i] for i in range(nthresh)]

	if labels is None:
		labels = ['p'+str(i) for i in range(nthresh)]

	if ax is None:
		fig, ax = plt.subplots(1	,1, constrained_layout=True)
	plt.grid()

	n, _, _ = create_hist_plot(data, nmin = minim, nmax = maxim, binsize = 1, ax = ax)


	for i in range(nthresh):

		mask = np.logical_and(domain > thresh_edges[i], domain < thresh_edges[i+1])
		subdomain = domain[mask]
		n_subset = n[mask]


		#fit damped sin
		popt, pcov = curve_fit(poissonian, subdomain, n_subset, p0 = [guesses[i],0.1])

		fit_values.append([popt, pcov])
		fitdomain = np.linspace(maxim, minim, 1000)

		ax.plot(fitdomain, poissonian(fitdomain, *popt), '--', c = colours[i], linewidth = 3, label = labels[i])

	ax.legend()
	ax.set_title('PMT counts')
	ax.set_xlabel('counts')
	ax.set_ylabel('frequency')

	return fit_values
