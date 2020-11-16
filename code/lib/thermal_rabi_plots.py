from lib import full_Rabi
import numpy as np 
import matplotlib.pyplot as plt

t = np.linspace(0,1000,1000)
#in microseconds  

#Now let out pi time be 30us (i.e. how well you can identify the carrier)
Rabi = np.pi/30
#Rabi now in MHz 

phase = 0 
#corresponds to axial frequency of 1MHz
dicke_factor = 0.068

n_av = [0.5,6,12,15,50]
det = [0, 0.5*Rabi, Rabi]

plt.rcParams.update({'font.size': 12})
fig_size = (5,3)
fig, axes = plt.subplots(fig_size[0]	,fig_size[1], constrained_layout=True)

#prepare figure



colours = ['b','m','c','r','tab:orange', 'tab:pink']


for i in range(fig_size[0]):
    for j in range(fig_size[1]):
        axes[i,j].plot(t, full_Rabi(t, [dicke_factor, n_av[i], det[j], Rabi, phase], m = 0), color = colours[i])
        axes[i,j].set_xticks(np.arange(0,1200,step = 200))
        axes[i,j].set_yticks(np.arange(0,1.2,step = 0.2))
        if j == 0:
            axes[i,j].set_ylabel("$n_{av} = $" + str(n_av[i]))
        if i == 0:
            axes[i,j].set_title("$\Delta = $" + str(round(det[j],3))+ "MHz")
fig.suptitle('Carrier', fontsize=16)
plt.savefig('C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/thermalrabicarrier.png', dpi = 1000)

plt.show() 

plt.rcParams.update({'font.size': 12})
fig_size = (5,3)
fig, axes = plt.subplots(fig_size[0]	,fig_size[1], constrained_layout=True)

for i in range(fig_size[0]):
    for j in range(fig_size[1]):
        axes[i,j].plot(t, full_Rabi(t, [dicke_factor, n_av[i], det[j], Rabi, phase], m = 1), color = colours[i])
        axes[i,j].set_xticks(np.arange(0,1200,step = 200))
        axes[i,j].set_yticks(np.arange(0,1.2,step = 0.2))
        if j == 0:
            axes[i,j].set_ylabel("$n_{av} = $" + str(n_av[i]))
        if i == 0:
            axes[i,j].set_title("$\Delta = $" + str(round(det[j],3))+ "MHz")
fig.suptitle('BSB', fontsize=16)
plt.savefig('C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/thermalrabiblue.png', dpi = 1000)

plt.show() 

plt.rcParams.update({'font.size': 12})
fig_size = (5,3)
fig, axes = plt.subplots(fig_size[0]	,fig_size[1], constrained_layout=True)

for i in range(fig_size[0]):
    for j in range(fig_size[1]):
        axes[i,j].plot(t, full_Rabi(t, [dicke_factor, n_av[i], det[j], Rabi, phase], m = -1), color = colours[i])
        axes[i,j].set_xticks(np.arange(0,1200,step = 200))
        axes[i,j].set_yticks(np.arange(0,1.2,step = 0.2))
        if j == 0:
            axes[i,j].set_ylabel("$n_{av} = $" + str(n_av[i]))
        if i == 0:
            axes[i,j].set_title("$\Delta = $" + str(round(det[j],3))+ "MHz")
fig.suptitle('RSB', fontsize=16)
plt.savefig('C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/thermalrabired.png', dpi = 1000)

plt.show() 
