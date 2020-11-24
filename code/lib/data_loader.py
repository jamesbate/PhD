"""This function is responsible for loading mulitple sets of data into
TricsDataObjects.
"""
##----------------------------PREAMBLE----------------------------------##
import numpy as np
import matplotlib.pyplot as plt
#external dependencies

from lib import TricsDataObject
#internal imports

##----------------------------MAIN----------------------------------##

def data_loader(filenames, data_folders, data_dir_trics, data_column = 6):
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
            data_objects[n,m] = trics_data_object

    return data_objects
