import numpy as np
import matplotlib.pyplot as plt
from lib import full_Rabi

n_av = 1
Rabi_0 = 0.934*2
dicke_factor = 0.06

#Let our ion be in some thermal state
# therm_prob = lambda n: (n_av/(n_av + 1))**n
# Rabi_dicke = lambda n: 0.5*(1 - dicke_factor**2*(n + 0.5))

t = np.linspace(670,710,10000)
#t = np.linspace(0,15,10000)


# def full_Rabi(t):
#     r = 0
#     for n in range(0,10):
#         r += therm_prob(n)*(np.sin(0.5*t*Rabi_dicke(n)))**2

#     return r
plt.title('Dicke Rabi Revivals!')
plt.plot(t, full_Rabi(t, [0.07, n_av, 0, Rabi_0, 0], m = 0))
plt.show()
