import numpy as np

def sense_qubit_rescale(s, nqubits_array):
    #for now just deal with one qubit in each sensor, this will shortly be extended
    s_max = np.amax(s,axis=0)
    s_scaled = s/s_max

    return s_scaled
