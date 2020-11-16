from lib import TricsDataObject, data_loader, create_hist_plot, rabi_fitter
import matplotlib.pyplot as plt
import numpy as np

##----------------------------PARAMETERS----------------------------------##
#data_dir_trics = "zidshare.uibk.ac.at/qos/qfc/measurements/trics_data/"
data_dir_trics = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/'
save_dir = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/trics/'

data_folders = ["2020-10-5/164939/"]
filenames = ["PMT1_2.txt"]


savefiles = ['testrab_exc_2020-10-5.png']

##----------------------------MAIN----------------------------------##

data_object = data_loader(filenames, data_folders, data_dir_trics)[0][0]
phases = data_object.df['seq.pulse729_length'].to_numpy()
(probs, dprobs) = data_object.get_probs(thresholds = [8])
initial_fits = [
[0.06       , 13.1, 0, 0.16632*np.pi - 0.19]
]

results = rabi_fitter(phases, np.array([probs[:,0]]).T, np.array([dprobs[:,0]]).T,labels = ['p2'], initial_fit = initial_fits, ground = True)

print('Rabi: ',results[0][3]/2)
print('Phonons: ',results[0][1])
print('lamb-dicke: ',results[0][0])
plt.savefig(save_dir + savefiles[0], dpi = 1000)
plt.show()
