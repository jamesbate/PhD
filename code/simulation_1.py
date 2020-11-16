from qutip import sigmax, sigmaz, mesolve, basis, sigmay
import matplotlib.pyplot as plt
import numpy as np

ham = 2*np.pi*0.1*sigmax()
psi_0 = basis(2,0)

times = np.linspace(0,10,20)

result = mesolve(ham,psi_0, times, [], [sigmaz(), sigmay()])


fig, ax = plt.subplots()
ax.plot(result.times, result.expect[0])
ax.plot(result.times, result.expect[1])

plt.show()
