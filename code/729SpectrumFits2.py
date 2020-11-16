from lib import TricsDataObject, data_loader, spectrum_lorentzian_fitter
import matplotlib.pyplot as plt
import numpy as np

##----------------------------PARAMETERS----------------------------------##
#data_dir_trics = "zidshare.uibk.ac.at/qos/qfc/measurements/trics_data/"
data_dir_trics = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/'
save_dir = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/trics/'

data_folders = ["2020-04-24/101858/"]
filenames = ["PMT1_2.txt"]

savefiles = ['729spectrumfit2_1.png']

##----------------------------MAIN----------------------------------##

[data_obj_1] = data_loader(filenames, data_folders, data_dir_trics)[:,0]
freq = data_obj_1.df['transitions.729_Carrier.shift'].to_numpy()
data_obj_1.histogram_plot(thresholds = [6])
plt.savefig(save_dir + "PMT"+savefiles[0], dpi = 1000)
plt.close()

(probs, dprobs) = data_obj_1.get_probs(thresholds = [6])
initial_p_carrier = [-0.3, 0.6, 0.05]
initial_p_radial = [1.8, 0.5, 0.05, -2.5, 0.5, 0.05]
params = spectrum_lorentzian_fitter(np.array([freq, freq]), np.array([probs[:,0], probs[:,0]]), np.array([dprobs[:,0], dprobs[:,0]]),initial_parameters =[initial_p_carrier, initial_p_radial], labels = ['carrier','radial'])

params = np.array(params)
print('Carrier:\n')
print('freq: ', params[0,0][1])
print('width:',params[0,0][3])

print('\nRadial Sidebands:\n')
print('freq: ', params[1,0][1])
print('width:',params[1,0][3])
print('dfreq: ', params[1,0][1] - params[0,0][1])
print('dfreqwidth: ', params[1,0][3] + params[0,0][3])
print('freq: ', params[1,0][4])
print('width:',params[1,0][6])
print('dfreq: ', params[1,0][4] - params[0,0][1])
print('dfreqwidth: ', params[1,0][6] + params[0,0][3])

print('\nRadial Average: ', np.mean([abs(params[1,0][1] - params[0,0][1]), abs(params[1,0][4] - params[0,0][1])]))
print(' pm ', np.mean([abs(params[1,0][3] + params[0,0][3]), abs(params[1,0][6] + params[0,0][3])]))

plt.savefig(save_dir +savefiles[0], dpi = 1000)
plt.close()
