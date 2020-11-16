import numpy as np
from qutip import Qobj
def tomo_matrix(t):
    return np.array([[t[0],0],[complex(t[2], t[3]), t[1]]])

def rho(t):
    tomo_mat = tomo_matrix(t)
    tomo_obj =  Qobj(tomo_mat)
    nrho = tomo_obj.dag()*tomo_obj/(tomo_obj.dag()*tomo_obj).tr()
    return nrho

print(rho([1,1,1,1]))
