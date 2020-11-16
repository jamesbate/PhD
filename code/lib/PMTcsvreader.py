"""Info
"""
##----------------------------PREAMBLE----------------------------------##
import numpy as np
from .trics_data_object_2 import TricsDataObject

##----------------------------PARAMETERS----------------------------------##
#data_dir_trics = "zidshare.uibk.ac.at/qos/qfc/measurements/trics_data/"
data_dir_trics = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/'

data_folders = ["171718/"]
filenames = ["PMT1_1.txt","PMT1_2.txt"]

##----------------------------MAIN----------------------------------##

def data_loader(filenames = filenames, data_folders = data_folders, data_dir_trics = data_dir_trics):
	"""Info
	"""
	#initialise
	folder_num = len(data_folders)
	file_num = len(filenames)
	data_objects = np.empty([folder_num,file_num], dtype = object)

	#column at which data starts
	data_column = 6

	for n,folder in enumerate(data_folders,0):
		for m,file in enumerate(filenames,0):

			#create object
			fields = np.genfromtxt(data_dir_trics + folder + file, delimiter="\t", max_rows = 1, dtype = str)
			trics_data_object = TricsDataObject(list(fields[:data_column]) + ['PMTcounts'], data_column)

			#load data
			data = np.genfromtxt(data_dir_trics + folder + file, delimiter="\t", skip_header = 1, dtype = float)
			trics_data_object.add_sequence_iter(data)

			data_objects[n][m] = trics_data_object

	return data_objects


if __name__ == "__main__":
	data_objects = data_loader()

	print(data_objects[0][0].column_names)
