import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

def create_hist_plot(data, nmin = 0, nmax = 100, binsize = 1 , thresholds = None, ax = None):

	binedges = np.arange(nmin - binsize/2, nmax + binsize/2, binsize)

	if ax is None:
		fig, ax = plt.subplots(1	,1, constrained_layout=True)
	n, bin_edges, _ = ax.hist(data, binedges, color = 'C2')

	if thresholds is not None:
		for t in thresholds:
			ax.axvline(x=t, color='r')

	if thresholds is not None:
		_fig = plt.figure()
		n, bin_edges, _ = plt.hist(data, [binedges[0]] + thresholds + [binedges[-1]])
		plt.close(_fig)


	return (n, n/data.size, bin_edges)
