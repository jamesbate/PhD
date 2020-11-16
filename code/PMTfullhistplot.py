from lib import data_loader, create_hist_plot, PMT_fitter
import matplotlib.pyplot as plt
import numpy as np

##----------------------------PARAMETERS----------------------------------##
#data_dir_trics = "zidshare.uibk.ac.at/qos/qfc/measurements/trics_data/"
data_dir_trics = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/'
save_dir = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/trics/'

data_folders = ["2020-10-5/164939/"]
filenames = ["PMT1_2.txt"]

savefiles = ['testhistogram2020-10-5.png']
##----------------------------MAIN----------------------------------##

data_objects = data_loader(filenames, data_folders, data_dir_trics)


for [data_object] in data_objects:

	[data] = data_object.get_sequences(['PMTcounts'])
	data_shaped = data.reshape(1,data_object.entries*len(data[0]))[0]
	create_hist_plot(data_shaped, nmin = 0, nmax = 100, thresholds = [8])

	result = PMT_fitter(data_shaped, thresholds = [8], labels = ['No ion', '1 ion'])
	plt.axvline(x = 8, c = 'r')
	print('fit parameters:\n', result)
	plt.savefig(save_dir + savefiles[0], dpi = 1000)
	plt.show()
