from .gram_schmidt import sub_gram_schmidt
import numpy as np

def sense_optimal_vectors(x):
    opt_vectors = np.zeros(x.shape)
    for i in range(x.shape[1]):
        opt_vectors[:,i] = sub_gram_schmidt(x[:,i], np.delete(x,i,axis = 1))
    return opt_vectors

if __name__ == "__main__":
    v = np.array([[1,1,1,1,1],[-2.12,-1,0,1,2.12],[4.4944,1,0,1,4.4944]])
    print(sense_optimal_vectors(v.T))
