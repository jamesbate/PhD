from lib import TricsDataObject, data_loader, rabi_fitter, rabi_fitter_badpump, plot_error_joint
import matplotlib.pyplot as plt
import numpy as np

##----------------------------PARAMETERS----------------------------------##
#data_dir_trics = "zidshare.uibk.ac.at/qos/qfc/measurements/trics_data/"
data_dir_trics = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/'
save_dir = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/trics/'

data_folders = ["2020-11-06/182939/", "2020-11-06/183239/","2020-11-06/180533/" ,"2020-11-06/181019/", "2020-11-06/180634/", "2020-11-06/180747/", "2020-11-06/181122/"]
data_folders_2 = ["2020-11-06/181343/","2020-11-06/182002/","2020-11-06/182100/","2020-11-06/182207/","2020-11-06/182316/"]
data_folder_3 = ['2020-11-06/183628/']
data_folder_4 = ['2020-11-06/182740/']
filenames = ["PMT1_2.txt"]
filenames_3 = ["camera_all.txt"]
 
savefiles = ['blue_temperature_fits2std.png', 'red_temperature_fits2std.png']

##----------------------------MAIN----------------------------------##
#prepare plot for joint plots
fig, ax = plt.subplots(1	,1, constrained_layout=True, figsize=(18, 9))
plt.grid()
colours = ['b','m','c','r','g','tab:orange', 'tab:pink']

bounds = [1.5,5.99,9]
i = 4



save = str(i)

#data_objs = data_loader(filenames, data_folders_2, data_dir_trics)[:,0]
data_objs = data_loader(filenames, data_folders, data_dir_trics)[2:,0]
do = data_objs[i]
pulse_length = do.df['seq.pulse729_length'].to_numpy() 
me, dme = do.get_mean_exc(thresholds = [18,75],  nmax = 150)
initial_fits = [
    [0.068,3,0,0.24375162841468767,0,0.9]
]
res = []
lab = []
for b in bounds:
    results = rabi_fitter_badpump(pulse_length, me, dme,labels = ['Mean Excitation'], initial_fit = initial_fits, p_fix =[[0.068,b,None,0.21*2,0,0.9]], m = 1,color = colours[4])
    res.append(results[0][1][2])
    lab.append(results)
    plt.close()

final_result = lab[1]
fit_label = "Rabi: " + str(round(final_result[0][3][0]/2,3)) + "MHz\nPhonons: " + str(round(final_result[0][1][0],3))+ ' ('+str(bounds[0])+ ', '+str(bounds[2])+') '+ "\nDicke factor: " + str(round(final_result[0][0][0],3))

plot_error_joint(pulse_length, me,dme,res[1],res[0], res[2],ax = ax, label = fit_label, color = colours[i])

plt.savefig(save_dir+save+savefiles[0], dpi = 1000)
#plt.savefig(save_dir+savefiles[1], dpi = 1000)
plt.show()