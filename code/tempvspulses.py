import numpy as np 
import matplotlib.pyplot as plt 
from lib import PhysicalConstants, phonon_to_temp, temp_to_phonon
save_dir = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/images/trics/'
pc = PhysicalConstants() 
freq = 1e6


plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots(1	,1, constrained_layout=True, figsize=(18, 9))

plt.grid()

#blue
raman_pulses = np.array([0,10,20,20,40])
phonon_upper = np.array([0.3,1.8,1.1,3,9])
phonons = np.array([0.14,0.72,1.3,1.56,5.99])
phonon_lower = np.array([0.01,0.6,0.4,0.9,1.5])
ax.scatter(raman_pulses, phonons, c = 'b')
ax.errorbar(raman_pulses, phonons, yerr = [phonon_lower, phonon_upper],ls = 'none', capsize = 3, c='b')


#red
raman_pulses = np.array([10,20,20,40,60])
phonon_upper = np.array([0.9,1.8,2,6,15])
phonons = np.array([0.67,1.15,1.22,3.22,9.41])
phonon_lower = np.array([0.4,0.7,0.8,2,5])
ax.scatter(raman_pulses, phonons, c = 'r')
ax.errorbar(raman_pulses, phonons, yerr = [phonon_lower, phonon_upper],ls = 'none', capsize = 3,c='r')



# secaxy = ax.secondary_yaxis('right', functions=(phonon_to_temp, temp_to_phonon))
# secaxy.set_ylabel(r'Temperature $(mK)$')

plt.title('Ion temperature vs Raman pulse number')
plt.xlabel('Pulses')
ax.set_ylabel(r'$\langle n\rangle$')                                                                      
ax2 = ax.twinx() 

ax2.set_ylabel('Temperature (mK)')
ax.set_ylim(0,25)
ax2.set_ylim(0,phonon_to_temp(np.array([25]), freq = 2*np.pi*1e6))
plt.savefig(save_dir+'tempvspulses.png', dpi = 1000)
plt.show() 
