"""This file contains functions useful for quantum sensing
"""
##----------------------------PREAMBLE----------------------------------##
import numpy as np
from lib import sub_gram_schmidt, gram_schmidt_columns

##------------------------FUCNTION DEFINITIONS------------------------------##

def sense_echo_times(s_scaled_array, nqubits_array):
    """calculate when to flip states given s array
    """
    #return fraction of sensing time
    echoes_array = np.empty(s_scaled_array.shape)
    initial_state_array = np.zeros(s_scaled_array.shape)

    i = 0
    for s_scaled, initial_state in zip(s_scaled_array.T, initial_state_array.T):
        echoes = echoes_array.T[i]
        echoes = 0.5*(1 - s_scaled/nqubits_array)
        for n in range(s_scaled.size):
            if echoes[n] > 0.5:
                echoes[n] = 1 - echoes[n]
                initial_state[n] = 1
            if echoes[n] < 10**-10:
                echoes[n] = 0
        echoes_array.T[i] = echoes
        i+=1
    return (echoes_array, initial_state_array)


def sense_optimal_vectors(x):
    """calculates s_opt using gram schmidt
    """
    opt_vectors = np.zeros(x.shape)
    for i in range(x.shape[1]):
        opt_vectors[:,i] = sub_gram_schmidt(x[:,i], np.delete(x,i,axis = 1))
    return opt_vectors
    

def sense_qubit_rescale(s, nqubits_array):
    """works out how to scale spin normalisation
    """
    #for now just deal with one qubit in each sensor, this will shortly be extended
    s_max = np.amax(s,axis=0)
    s_scaled = s/s_max

    return s_scaled
