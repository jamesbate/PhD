import numpy as np
from py_pol.jones_vector import Jones_vector, degrees
from py_pol.stokes import Stokes

S = Stokes('Linear light')
S.linear_light(azimuth=np.linspace(0,180,13)*degrees)
ax, fig = S.draw_poincare()
