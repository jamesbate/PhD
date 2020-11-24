"""This file contains functions which do full and sub gram schmidt routines
"""
##-------------------------------PREAMBLE-----------------------------------##
import numpy as np
#external dependencies

##-------------------------------PARAMETERS-----------------------------------##

#below this value, floats are taken to be zero if zero_threshold flag is true
zero_threshold = 1e-10

##-------------------------FUNCTION DEFINITIONS-----------------------------##

def gram_schmidt_columns(X, zero_threshold = True):
    """Full gram schmidt procedure on columns of matrix
    Parameters
    ----------
    X : numpy matrix
        matrix of column vectors

    zero_threshold : bool
        below treshold, all floats taken to be zero

    Returns
    -------
    numpy matrix
        matrix after gram schmidt
    """
    Q, R = np.linalg.qr(X)

    if zero_threshold:
        #arbitrarily set 0 threshold
        mask = abs(Q) < zero_threshold
        Q[mask] = 0

    return Q

def sub_gram_schmidt(X_to_adjust, X_rest, zero_threshold = True):
    """Single step in gram schmidt procedure
    Parameters
    ----------
    X_to_adjust : numpy array
        column vector to make normal to rest

    X_rest : numpy matrix
        matrix of rest of column vectors

    zero_threshold : bool
        below treshold, all floats taken to be zero

    Returns
    -------
        vector which is now normal to rest
    """
    #noramlise everything
    X_rest_sums = sum(X_rest**2)
    X_rest = X_rest / np.sqrt(X_rest_sums)
    X_to_adjust = X_to_adjust/np.sqrt(sum(X_to_adjust**2))

    for v in X_rest.T:
        X_to_adjust -= np.dot(X_to_adjust, v)*v

    if zero_threshold:
        #arbitrarily set 0 threshold
        mask = abs(X_to_adjust) < zero_threshold
        X_to_adjust[mask] = 0

    return X_to_adjust

if __name__ == "__main__":
    v = np.array([[1,1,1,1,1],[-2.12,-1,0,1,2.12],[4.4944,1,0,1,4.4944]])
    #tranmspose as must use column vectors 
    print(sub_gram_schmidt(v[0].T, v[1:3].T))
