"""Creates histogram plot of data with verticle threshold markers. If thresholds
are given, also returns populations within each bin.
"""
##----------------------------PREAMBLE----------------------------------##
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
#external dependencies

##------------------------FUNCTION DELCARATION------------------------------##

def plot_histogram(data, nmin = 0, nmax = 100, binsize = 1 , thresholds = None, ax = None):
	"""
	Parameters
	----------
	data : array-like
		data for histogram

	nmin : int
		smallest bin edge

	nmax : int
		largest bin edge

	binsize : int
		sizes of bins

	thresholds : array-like
		thresholds to distinguish different numbers of ions

	ax : mpl figure
		plt onto another axis

	Returns
	-------
	(n, n/data.size, bin_edges) : (array-like, int, array-like)
		n is the populations of each bin, second entry is probability after
		normalising 
	"""
	binedges = np.arange(nmin - binsize/2, nmax + binsize/2, binsize)

	#prepare figure
	if ax is None:
		fig, ax = plt.subplots(1	,1, constrained_layout=True)
	n, bin_edges, _ = ax.hist(data, binedges, color = 'C2')

	#plot verticle lines for thresholds
	if thresholds is not None:
		for t in thresholds:
			ax.axvline(x=t, color='r')

		#temporarilty use plt.hist to bin data. If binning becomes more common
		#will have to create binning function in lib
		_fig = plt.figure()
		n, bin_edges, _ = plt.hist(data, [binedges[0]] + thresholds + [binedges[-1]])
		plt.close(_fig)

	return (n, n/data.size, bin_edges)
