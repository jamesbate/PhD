import matplotlib.pyplot as plt 
import numpy as np 
from lib import dicke_regime_flops_blue, dicke_regime_flops_red, dicke_regime_flops_carrier, full_Rabi, full_rabi_dicke

#prepare figure
plt.rcParams.update({'font.size': 16})
colours = ['b','m','c','r','tab:orange', 'tab:pink']
fig, axes = plt.subplots(1	,2, constrained_layout=True, figsize=(18, 9))
plt.grid()

t = np.linspace(0,1000,1000)
#in microseconds  

#Now let out pi time be 30us (i.e. how well you can identify the carrier)
Rabi = np.pi/50
#Rabi now in MHz 

#corresponds to axial frequency of 1MHz
eta = 0.068

n_av = 20

#These are the plots in the lambe dicke regime (up to 2nd order)
axes[0].plot(t, dicke_regime_flops_blue(t, Rabi, eta, n_av), ':', c = colours[0])
axes[0].plot(t, dicke_regime_flops_red(t, Rabi, eta, n_av), ':', c = colours[1])
axes[0].plot(t, dicke_regime_flops_carrier(t, Rabi, eta, n_av), ':', c = colours[2])

#This is the full solution
axes[0].plot(full_Rabi(t, [eta, n_av, 0, Rabi, 0] , m = 1), c = colours[0], label = "BSB")
axes[0].plot(full_Rabi(t, [eta, n_av, 0, Rabi, 0] , m = 0), c = colours[2], label = "Carrier")
axes[0].plot(full_Rabi(t, [eta, n_av, 0, Rabi, 0] , m = -1), c = colours[1], label = "RSB")

axes[0].set_title("$n_{av} = 20$")

n_av = 0.5

#These are the plots in the lambe dicke regime (up to 2nd order)
axes[1].plot(t, dicke_regime_flops_blue(t, Rabi, eta, n_av), ':', c = colours[0])
axes[1].plot(t, dicke_regime_flops_red(t, Rabi, eta, n_av), ':', c = colours[1])
axes[1].plot(t, dicke_regime_flops_carrier(t, Rabi, eta, n_av), ':', c = colours[2])

#This is the full solution
axes[1].plot(full_Rabi(t, [eta, n_av, 0, Rabi, 0] , m = 1), c = colours[0], label = "BSB")
axes[1].plot(full_Rabi(t, [eta, n_av, 0, Rabi, 0] , m = 0), c = colours[2], label = "Carrier")
axes[1].plot(full_Rabi(t, [eta, n_av, 0, Rabi, 0] , m = -1), c = colours[1], label = "RSB")

axes[1].set_title("$n_{av} = 0.5$")

plt.legend()
plt.savefig('C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/dickevsfabiflops.png', dpi = 1000)
plt.show()