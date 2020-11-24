"""This file contains functions useful for state tomography  
"""
##---------------------------PREAMBLE-----------------------------------------##
from qutip import expect, sigmax, sigmaz, sigmay, qeye, tensor
import numpy as np
##----------------------------STOKES PARAMETERS----------------------------##

def generalised_stokes_projection(rho, stokes_indices, bases = [qeye(2), sigmax(), sigmay(), sigmaz()]):
    exp = []
    dim = stokes_indices.shape[1]
    for indices in stokes_indices:
        exp.append(expect(rho, tensor([bases[index] for index in indices])))
    return exp

def bases_rotation_matrix():
    pass
