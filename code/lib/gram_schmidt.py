import numpy as np

def gram_schmidt_columns(X, zero_threshold = True):
    Q, R = np.linalg.qr(X)

    if zero_threshold:
        #arbitrarily set 0 threshold at 10^-10
        mask = abs(Q) < 1e-10
        Q[mask] = 0
    return Q

def sub_gram_schmidt(X_to_adjust, X_rest, zero_threshold = True):
    #noramlise everything
    X_rest_sums = sum(X_rest**2)
    X_rest = X_rest / np.sqrt(X_rest_sums)
    X_to_adjust = X_to_adjust/np.sqrt(sum(X_to_adjust**2))

    for v in X_rest.T:
        X_to_adjust -= np.dot(X_to_adjust, v)*v

    if zero_threshold:
        #arbitrarily set 0 threshold at 10^-10
        mask = abs(X_to_adjust) < 1e-10
        X_to_adjust[mask] = 0

    return X_to_adjust

if __name__ == "__main__":
    #give the function the example from sensing paper
    #x = np.array([[1,1,1,1,1],[4,1,0,1,4],[64,1,0,1,64],[-2,-1,0,1,2],[-8,-1,0,1,8]]).T

    v = np.array([[1,1,1,1,1],[-2.12,-1,0,1,2.12],[4.4944,1,0,1,4.4944]])
    print(sub_gram_schmidt(v[0].T, v[1:3].T))
