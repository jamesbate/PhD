##---------------------------PREAMBLE-----------------------------------------##
from qutip import ket2dm, basis, Qobj, tensor
import numpy as np
from lib import generalised_stokes_projection
##----------------------------STOKES PARAMETERS----------------------------##

rho_H = ket2dm(basis(2,0))
rho = Qobj([[5/8,complex(0,-1/(2*np.sqrt(2)))],[complex(0,1/(2*np.sqrt(2))), 3/8]])
rho_2d = tensor(ket2dm(basis(2,0)), ket2dm(basis(2,0)))

print(generalised_stokes_projection(rho_2d, np.array([[0,0],[3,3], [0,1]])))
print(generalised_stokes_projection(rho, np.array([[0],[1],[2]])))
