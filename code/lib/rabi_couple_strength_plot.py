from lib import full_rabi_dicke
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.special import eval_genlaguerre, factorial

n = np.arange(0,10,1)

plt.plot(n, [full_rabi_dicke(i,0, 0.05, 1) for i in n], label = "Carrier")
plt.plot(n, [full_rabi_dicke(i,-1, 0.05, 1) for i in n], label = "Red Sideband")
plt.plot(n, [full_rabi_dicke(i,1, 0.05, 1) for i in n], label = "Blue Sideband")
plt.plot(n, [full_rabi_dicke(i,2, 0.05, 1) for i in n], label = "2nd Sideband")
plt.plot(n, [full_rabi_dicke(i,3, 0.05, 1) for i in n], label = "3rd Sideband")

plt.plot(n, 1*0.05*np.sqrt(n + 1), ':', label = "1st Sideband LD")

plt.ylabel("Relative coupling strength")
plt.xlabel("n")
plt.ylim(top = 1)
plt.legend()
plt.savefig('C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/rabidickefullclose.png', dpi = 1000)
plt.show() 