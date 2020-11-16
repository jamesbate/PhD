"""Info
"""
import numpy as np
import matplotlib.pyplot as plt
from .create_hist_plot import create_hist_plot

class TricsDataObject:
	"""Info
	data in form of dictionary of lists
	"""
	def __init__(self, column_names, data_column):
		self.column_names = column_names
		self.data = {}
		self._initialise_data_dict()
		self.entries = 0
		self.data_column = data_column
		self.data_name = column_names[data_column]

	def _initialise_data_dict(self):
		for c in self.column_names:
			self.data.update({c: []})


	def add_sequence(self, new_data_dict):
		keys = new_data_dict.keys()
		if keys != self.data.keys():
			raise ValueError('data mismatch, enter values for all entries!')
		for k in keys:
			self.data[k] += [new_data_dict[k]]
		self.entries += 1

	def add_sequence_list(self, data_list):

		for n in range(self.data_column):
			self.data[self.column_names[n]] += [data_list[n]]

		self.data[self.data_name] += [data_list[self.data_column:]]
		self.entries += 1

	def add_sequence_iter(self, data):
		for d in data:
			self.add_sequence_list(d)

	def get_sequences(self, column_names = None, indices = None):
		if column_names is None:
			column_names = self.column_names
		if indices is None:
			indices = list(range(self.entries))

		data_subset = []
		for col in column_names:
			d = [self.data[col][i] for i in indices]
			data_subset += [d]
		return np.array(data_subset)

	def sample_histograms_sim(self, indices = list(range(5)), thresholds = None):
		#assumed from self.data_name

		nplots = len(indices)

		#extract data
		[data] = self.get_sequences([self.data_name], indices)

		#make hist plots
		fig, axes = plt.subplots(nplots	,1, constrained_layout=True)

		for m,i in enumerate(indices,0):
			n, _ = create_hist_plot(data[i], np.arange(0,100,1), ax = axes[m])

		if thresholds is not None:
			for t in thresholds:
				for ax in axes:
					ax.axvline(x=t, c = 'r')

	def full_histogram_plot(self, thresholds = None):
		#assumed from self.data_name

		#extract data
		[data] = self.get_sequences([self.data_name])
		data_shaped = data.reshape(1,self.entries*len(data[0]))[0]

		#make hist plots
		fig, ax = plt.subplots(1	,1, constrained_layout=True)
		plt.grid()
		_, p, _ = create_hist_plot(data_shaped, nmin = 0, nmax = 100, binsize = 1, ax = ax)

		if thresholds is not None:
			for t in thresholds:
				ax.axvline(x=t, c = 'r')
		return n

	def get_probs(self, thresholds, indices = None):
		[data] = self.get_sequences([self.data_name], indices)
		n = np.zeros([self.entries,len(thresholds) + 1])

		for i,d in enumerate(data,0):
			_, _p, _ = create_hist_plot(d, nmin = 0, nmax = 100, binsize = 1, thresholds = thresholds)
			n[i] = _p
			plt.close()

		return n
