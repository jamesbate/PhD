import numpy as np

def sense_echo_times(s_scaled_array, nqubits_array):
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
