from lib import TricsDataObject, data_loader, rabi_fitter, rabi_fitter_badpump
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
 
savefiles = ['blue_temperature_fits2_rabiset.png', 'red_temperature_fits2_rabiset.png']


##----------------------------MAIN----------------------------------##

#prepare plot for joint plots
fig, ax = plt.subplots(1	,1, constrained_layout=True, figsize=(18, 9))
plt.grid()
colours = ['b','m','c','r','g','tab:orange', 'tab:pink']

#load in data
data_objs = data_loader(filenames, data_folders, data_dir_trics)[:,0]


#FINAL CARRIER
# data_again = data_loader(filenames, data_folder_4, data_dir_trics)[0,0]
# pulse_length = data_again.df['seq.pulse729_length'].to_numpy()
# (probs, dprobs) = data_again.get_probs(thresholds = [18,70])
# me, dme = data_again.get_mean_exc(thresholds = [18,70], nmax = 150)

# #[dicke_factor, n_av, detuning, Rabi, phase, amp]
# initial_fits = [
# [0.068       , 1, 0, 0.31, 1, 0.9]
# ]
# results = rabi_fitter_badpump(pulse_length, np.array(me), np.array(dme),labels = ['Mean Excitation'], initial_fit = initial_fits, p_fix =[[0.068,None,None,None,None,None]], m = 0, ax = ax)
# print(results)
# #plt.savefig(save_dir + "beforerepump"+savefiles[0], dpi = 1000)




#CAMERA
camera_data = data_loader(filenames_3, data_folder_3, data_column=1)[0][0]
data = (camera_data.get_data(reshape=False))
scan_point = camera_data.df['# scan point'].to_numpy()
initial_fits = [
[0.068       , 0.1, 0, 0.47, 1],
[0.068       , 0.1, 0, 0.47, 1]
]
results = rabi_fitter(scan_point, np.array(data), np.array([[None], [None]]).T,labels = ['p_1', 'p_1'], initial_fit = initial_fits, p_fix =[[0.068,None,0,None,0.83]]*2, m = 0, ax = ax)
print(results)
plt.savefig(save_dir + "camera"+savefiles[0], dpi = 1000)
plt.show()
exit()



#by viewing histogram plot, I extract thresholds of [18,70]. I use nmax = 150.
#data_obj_1.histogram_plot(thresholds = [18, 70],nmax = 150) 


# pulse_length = data_objs[1].df['seq.pulse729_length'].to_numpy() 
# (probs, dprobs) = data_objs[1].get_probs(thresholds = [18,70])
# me, dme = data_objs[1].get_mean_exc(thresholds = [18,70], nmax= 150)
# #[dicke_factor, n_av, detuning, Rabi, phase]
# initial_fits = [
# [0.068       , 6, 0, 0.31, 1]
# ]
# results = rabi_fitter(pulse_length, np.array(me), np.array(dme),labels = ['Mean Excitation'], initial_fit = initial_fits, p_fix =[[0.068,None,None,None,None]], m = 0, ax = ax, color = ['g'])
# print(results)
# #plt.savefig(save_dir + "afterrepump"+savefiles[0], dpi = 1000)

# plt.close()




# pulse_length = data_objs[0].df['seq.pulse729_length'].to_numpy()
# (probs, dprobs) = data_objs[0].get_probs(thresholds = [18,70])
# me, dme = data_objs[0].get_mean_exc(thresholds = [18,70], nmax = 150)

# #[dicke_factor, n_av, detuning, Rabi, phase, amp]
# initial_fits = [
# [0.068       , 6, 0, 0.31, 1, 0.9]
# ]
# results = rabi_fitter_badpump(pulse_length, np.array(me), np.array(dme),labels = ['Mean Excitation'], initial_fit = initial_fits, p_fix =[[0.068,None,None,None,None,None]], m = 0, ax = ax, color = ['r'])
# print(results)
# #plt.savefig(save_dir + "beforerepump"+savefiles[0], dpi = 1000)
# plt.show()
# plt.close() 
# exit()




# for i in range(5):
#     do = data_objs[i+2]

#     pulse_length = do.df['seq.pulse729_length'].to_numpy() 
#     me, dme = do.get_mean_exc(thresholds = [18,75], nmax = 150)
#     initial_fits = [
#     [0.068       , 2, 0, np.pi/8, 0, 0.9]
#     ]
#     results = rabi_fitter_badpump(pulse_length, me, dme,labels = ['Mean Excitation'], initial_fit = initial_fits, p_fix =[[0.068,None,0,2*0.24,0,0.9]], m = 1,color = colours[i])
#     print(results)

#     plt.savefig(save_dir + str(i)+savefiles[0], dpi = 1000)
# # plt.savefig(save_dir + savefiles[0], dpi = 800)
# # plt.show()
# # exit()

# #

# #Now look at red sidebands
# #load in data
# data_objs = data_loader(filenames, data_folders_2, data_dir_trics)[:,0]

# for i in range(3):
#     do = data_objs[i]
#     pulse_length = do.df['seq.pulse729_length'].to_numpy() 
#     me, dme = do.get_mean_exc(thresholds = [18,75],  nmax = 150)
#     initial_fits = [
#     [0.068       , 0.5, 0, np.pi/8, 0, 0.9]
#     ]
#     results = rabi_fitter_badpump(pulse_length, me, dme,labels = ['Mean Excitation'], initial_fit = initial_fits, p_fix =[[0.068,None,0,2*0.24,None,0.9]], m = -1,color = colours[i])


#     plt.savefig(save_dir + str(i)+savefiles[1], dpi = 1000)
#     print(results)


# do = data_objs[3]
# pulse_length = do.df['seq.pulse729_length'].to_numpy() 
# me, dme = do.get_mean_exc(thresholds = [18,75],  nmax = 150)
# initial_fits = [
# [0.068       , 5, 0, np.pi/8, 0, 0.9]
# ]
# results = rabi_fitter_badpump(pulse_length, me, dme,labels = ['Mean Excitation'], initial_fit = initial_fits, p_fix =[[0.068,None,0,2*0.24,0,0.9]], m = -1,color = colours[3])
# plt.savefig(save_dir+ str(3)+savefiles[1], dpi = 1000)
# print(results)

# do = data_objs[4]
# pulse_length = do.df['seq.pulse729_length'].to_numpy() 
# me, dme = do.get_mean_exc(thresholds = [18,75],  nmax = 150)
# initial_fits = [
# [0.068       , 7, 0, 0.7*np.pi/8, 0, 0.9]
# ]
# results = rabi_fitter_badpump(pulse_length, me, dme,labels = ['Mean Excitation'], initial_fit = initial_fits, p_fix =[[0.068,None,0,2*0.24,0,0.9]], m = -1,color = colours[4])
# plt.savefig(save_dir+ str(4)+savefiles[1], dpi = 1000)
# #plt.savefig(save_dir+savefiles[1], dpi = 1000)
# plt.show()
# print(results)

