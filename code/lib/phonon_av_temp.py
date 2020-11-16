import numpy as np 
from lib import PhysicalConstants

pc = PhysicalConstants() 


def phonon_to_temp(phonons, freq = 1e6, mk = True):
    l = phonons.size
    if np.any(np.logical_not(phonons > 0)):
        return [0]*l
    if l == 0:
        return []
    return 1e3*pc.h_bar*freq/(pc.k_B*np.log(1+1/phonons))

def temp_to_phonon(temp, freq = 1e6, mk = True):
    pass