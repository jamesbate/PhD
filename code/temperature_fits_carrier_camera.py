from lib import TricsDataObject, data_loader, rabi_fitter, rabi_fitter_badpump
import matplotlib.pyplot as plt
import numpy as np

##----------------------------PARAMETERS----------------------------------##
#data_dir_trics = "zidshare.uibk.ac.at/qos/qfc/measurements/trics_data/"
data_dir_trics = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/'
save_dir = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/trics/'

data_folders = ["2020-11-06/183628/"]
#data_folders_2 = ["2020-11-06/181343/","2020-11-06/182002/","2020-11-06/182100/","2020-11-06/182207/","2020-11-06/182316/"]
filenames = ["camera_ion1.txt", "camera_ion2.txt"]

savefiles = ['carriertempcam1.png', 'carriertempcam2.png']


##----------------------------MAIN----------------------------------##

#prepare plot for joint plots
fig, ax = plt.subplots(1	,1, constrained_layout=True, figsize=(18, 9))
plt.grid()
colours = ['b','m','c','r','tab:orange', 'tab:pink']

#load in data
data_objs = data_loader(filenames, data_folders, data_dir_trics, 3)[:,0][0]

scan = data_objs.df['# scan point']
print(data_objs.get_data(reshape = False))
