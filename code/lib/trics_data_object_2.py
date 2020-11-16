import numpy as np
import matplotlib.pyplot as plt
from .create_hist_plot import create_hist_plot
import pandas as pd
from collections import Counter

class TricsDataObject:
    def __init__(self, filename, data_start_column):
        self.df = self.load_data(filename)
        self.data_start_column = data_start_column
        self.column_names = self.df.columns

    def load_data(self, filename):
        return pd.read_csv(filename, delimiter = "\t")

    def get_data(self, reshape = True, indices = None):
        if indices is not None:
            data = self.df.iloc[indices, self.data_start_column:self.df.shape[1]].to_numpy()
        else:
            data = self.df.iloc[:,self.data_start_column:].to_numpy()
        if reshape == False:
            #transpose
            return data
        reshaped = np.prod(data.shape)
        return data.reshape(1,reshaped)[0]

    def histogram_plot(self, thresholds = None, indices = None, nmin = 0, nmax = 100):
    	#assumed from self.data_name

    	#extract data
        if indices is not None:
            mask = self.df.index.isin(indices)
            self.df = self.df[mask]

        data_shaped = self.get_data()

        #make hist plots
        fig, ax = plt.subplots(1	,1, constrained_layout=True)
        plt.grid()
        _, n, _ = create_hist_plot(data_shaped, nmin = nmin, nmax = nmax, binsize = 1, ax = ax)

        if thresholds is not None:
        	for t in thresholds:
        		ax.axvline(x=t, c = 'r')
        return n

    def get_probs(self, thresholds, indices = None, nmin = 0, nmax = 100):

        if indices is not None:
            data = self.get_data(reshape = False, indices = indices)
        else:
            data = self.get_data(reshape = False)

        p = np.zeros([data.shape[0],len(thresholds) + 1])
        dp = np.zeros([data.shape[0],len(thresholds) + 1])

        for i,d in enumerate(data,0):
            n = self.bin_data(d, nmin = nmin, nmax = nmax, binsize = 1, thresholds = thresholds)
            p[i] = np.array(n)/sum(n)
            dp[i] = np.maximum(np.sqrt(p[i]*(1-p[i])/sum(n)), 1/(sum(n)+2))
        return (p,dp)
    
    def get_mean_exc(self, thresholds, indices = None, nmin = 0, nmax = 100):
        (probs, dprobs) = self.get_probs(thresholds, indices = None, nmin = 0, nmax = 100)
        nions = len(thresholds)
        mean_exc = np.zeros([probs.shape[0],1])
        dprob_mean= np.zeros([probs.shape[0],1])
        for n in range(nions + 1):
            p = probs[:,n].reshape((probs.shape[0],1))
            dp = dprobs[:,n].reshape((probs.shape[0],1))
            mean_exc += n*p/nions
            dprob_mean += n*dp/nions

        return (1 - mean_exc, dprob_mean)


    def bin_data(self, data, thresholds, nmin = 0, nmax = 100, binsize = 1):
        cnt = Counter()
        bin_edges = np.array([nmin - binsize/2] + list(thresholds) + [nmax + binsize/2])
        bins = bin_edges[:-1]
        #initialise counter
        for e in bins:
            cnt[e] = 0
        #bin count data
        for d in data:
            for j in range(len(bin_edges)):
                if bin_edges[j] > d:
                    cnt[bins[j-1]] += 1
                    break
        return list(cnt.values())

if __name__ == "__main__":
    test_trics = TricsDataObject('C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/171718/PMT1_2.txt', 6)
    n = test_trics.get_probs([15])
    print(n)
