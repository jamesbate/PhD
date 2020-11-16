"""Info
"""
##----------------------------PREAMBLE----------------------------------##
import numpy as np
from . import TricsDataObject

##----------------------------PARAMETERS----------------------------------##
#data_dir_trics = "zidshare.uibk.ac.at/qos/qfc/measurements/trics_data/"
data_dir_trics = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/'

data_folders = ["171718/"]
filenames = ["PMT1_1.txt","PMT1_2.txt"]

##----------------------------MAIN----------------------------------##

def data_loader(filenames = filenames, data_folders = data_folders, data_dir_trics = data_dir_trics, data_column = 6):
    """The point of data_loader is to return an array of trics object for different sets of data, rather
    than doing each one manually.
    """

    folder_num = len(data_folders)
    file_num = len(filenames)
    data_objects = np.empty([folder_num,file_num], dtype = object)


    for n,folder in enumerate(data_folders,0):
    	for m,file in enumerate(filenames,0):

    		#create object
            trics_data_object = TricsDataObject(data_dir_trics + folder + file, data_column)
            data_objects[n][m] = trics_data_object

    return data_objects

if __name__ == "__main__":
    objects_array = data_loader()
    for object in objects_array:
        object.histogram_plot()
        plt.show()
