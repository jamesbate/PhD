from .gram_schmidt import gram_schmidt_columns
import numpy as np

def sense_vector(noise_signal_functions, sense_signal_functions, r):

    #project onto discrete spatial positions
    noise_signal_vectors = np.array([f(r) for f in noise_signal_functions]).T
    sense_signal_vectors = np.array([f(r) for f in sense_signal_functions]).T

    #this ordering is important! Must do noise before signal
    signal_vectors = np.concatenate((noise_signal_vectors, sense_signal_vectors), axis = 1)
    print('B:\n',signal_vectors)
    signal_vectors_gs = gram_schmidt_columns(signal_vectors)
    print('A:\n', signal_vectors_gs)

    return signal_vectors_gs[:,-1]

    #Now think about how to scale for available qubits

if __name__ == "__main__":

    #Taylor Generator
    @np.vectorize
    def taylor_generator(n):
        return lambda r: r**n

    #positions
    r = np.array([-2,-1,0,1,2])

    #noise signal functions matrix
    noise_signal_functions = taylor_generator([0,1,2,4])
    sense_signal_functions = taylor_generator([3])

    print(sense_vector(noise_signal_functions, sense_signal_functions, r))
