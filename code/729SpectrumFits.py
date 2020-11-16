from lib import TricsDataObject, data_loader, spectrum_lorentzian_fitter
import matplotlib.pyplot as plt
import numpy as np

##----------------------------PARAMETERS----------------------------------##
#data_dir_trics = "zidshare.uibk.ac.at/qos/qfc/measurements/trics_data/"
data_dir_trics = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/'
save_dir = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/trics/'

data_folders = ["2020-10-5/170354/", "2020-10-5/171540/"]
filenames = ["PMT1_2.txt"]

savefiles = ['729spectrumfit_1.png', '729spectrumfit_2.png']

##----------------------------MAIN----------------------------------##

[data_obj_1, data_obj_2] = data_loader(filenames, data_folders, data_dir_trics)[:,0]

freq = data_obj_1.df['transitions.729_Probe.frequency'].to_numpy()
data_obj_1.histogram_plot(thresholds = [6])
plt.savefig(save_dir + "PMT"+savefiles[0], dpi = 1000)
plt.close()

(probs, dprobs) = data_obj_1.get_probs(thresholds = [6])
inital_p = [451.9, 0.6, 0.1,450.9, 0.4, 0.1, 449.6, 0.5, 0.05,452.7, 0.5, 0.1, 453.9, 0.5, 0.05, 449.9, 0.5, 0.1]
initial_p_carrier = [451.9, 0.6, 0.1]
initial_p_axial = [450.9, 0.4, 0.1, 452.7, 0.5, 0.1, 449.9, 0.5, 0.1 ]
initial_p_radial = [453.9, 0.5, 0.05, 449.6, 0.5, 0.05]
params = spectrum_lorentzian_fitter(np.array([freq, freq, freq]), np.array([probs[:,0], probs[:,0], probs[:,0]]), np.array([dprobs[:,0], dprobs[:,0], dprobs[:,0]]),initial_parameters =[initial_p_carrier, initial_p_axial,initial_p_radial], labels = ['carrier','axial', 'radial'])

params = np.abs(params)
print('Carrier:\n')
print('freq: ', params[0,0][1])
print('width:',params[0,0][3])

print('\nAxial Sidebands:\n')
print('freq: ', params[1,0][1])
print('width:',params[1,0][3])
print('dfreq: ', params[1,0][1] - params[0,0][1])
print('dfreqwidth: ', params[1,0][3] + params[0,0][3])
print('freq: ', params[1,0][4])
print('width:',params[1,0][6])
print('dfreq: ', params[1,0][4] - params[0,0][1])
print('dfreqwidth: ', params[1,0][6] + params[0,0][3])
print('freq: ', params[1,0][7])
print('width:',params[1,0][9])
print('dfreq: ', params[1,0][7] - params[0,0][1])
print('dfreqwidth: ', params[1,0][9] + params[0,0][3])

print('\nAxial Average: ', np.mean([abs(params[1,0][1] - params[0,0][1]), abs(params[1,0][4] - params[0,0][1]), 0.5*(abs(params[1,0][7] - params[0,0][1]))]))
print(' pm ', np.mean([abs(params[1,0][3] + params[0,0][3]), abs(params[1,0][6] + params[0,0][3]), 0,5*(abs(params[1,0][9] + params[0,0][3]))]))

print('\nRadial Sidebands:\n')
print('freq: ', params[2,0][1])
print('width:',params[2,0][3])
print('dfreq: ', params[2,0][1] - params[0,0][1])
print('dfreqwidth: ', params[2,0][3] + params[0,0][3])
print('freq: ', params[2,0][4])
print('width:',params[2,0][6])
print('dfreq: ', params[2,0][4] - params[0,0][1])
print('dfreqwidth: ', params[2,0][6] + params[0,0][3])


print('n\Radial Average: ', np.mean([abs(params[2,0][1] - params[0,0][1]), abs(params[2,0][4] - params[0,0][1])]))
print(' pm ', np.mean([abs(params[2,0][3] + params[0,0][3]), abs(params[2,0][6] + params[0,0][3])]))


plt.savefig(save_dir + savefiles[0], dpi = 1000)
plt.close()

freq = data_obj_2.df['transitions.729_Probe.frequency'].to_numpy()
data_obj_2.histogram_plot(thresholds = [6])
plt.savefig(save_dir + "PMT"+savefiles[1], dpi = 1000)
plt.close()


(probs, dprobs) = data_obj_2.get_probs(thresholds = [6])
inital_p = [451.9, 0.6, 0.1,450.9, 0.4, 0.1, 449.6, 0.5, 0.05,452.7, 0.5, 0.1, 453.9, 0.5, 0.05, 449.9, 0.5, 0.1]
initial_p_carrier = [451.7, 0.6, 0.1]
initial_p_axial = [451.2, 0.4, 0.1, 452.3, 0.5, 0.1]
#initial_p_radial = [453.9, 0.5, 0.05, 449.6, 0.5, 0.05]
params = spectrum_lorentzian_fitter(np.array([freq, freq]), np.array([probs[:,0], probs[:,0]]), np.array([dprobs[:,0], dprobs[:,0]]),initial_parameters =[initial_p_carrier, initial_p_axial], labels = ['carrier','axial'])

plt.savefig(save_dir + savefiles[1], dpi = 1000)
plt.close()
params = np.abs(params)
print('\nCarrier:\n')
print('freq: ', params[0,0][1])
print('width:',params[0,0][3])

print('\nAxial Sidebands:\n')
print('freq: ', params[1,0][1])
print('width:',params[1,0][3])
print('dfreq: ', params[1,0][1] - params[0,0][1])
print('dfreqwidth: ', params[1,0][3] + params[0,0][3])
print('freq: ', params[1,0][4])
print('width:',params[1,0][6])
print('dfreq: ', params[1,0][4] - params[0,0][1])
print('dfreqwidth: ', params[1,0][6] + params[0,0][3])

print('\nAxial Average: ', np.mean([abs(params[1,0][1] - params[0,0][1]), abs(params[1,0][4] - params[0,0][1])]))
print(' pm ', np.mean([abs(params[1,0][3] + params[0,0][3]), abs(params[1,0][6] + params[0,0][3])]))
