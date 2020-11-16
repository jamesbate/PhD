##------------------preamble---------------------aa
import numpy as np 
import matplotlib.pyplot as plt 

##------------------parameters-------------------##

g_0 = lambda gamma,phi: 0.5*abs(np.cos(gamma)*np.sin(2*phi))
@np.vectorize 
def g_1(gamma,phi):
	return abs(complex(np.cos(gamma)*np.cos(2*phi), np.sin(gamma)*np.cos(phi)))/np.sqrt(6)  

@np.vectorize
def g_2(gamma,phi):
	return abs(complex(0.5*np.cos(gamma)*np.sin(2*phi), np.sin(gamma)*np.sin(phi)))/np.sqrt(6)

##-----------------plot---------------------------##

#The grid of the plot must be manually adjusted

if __name__ == "__main__":
	#create meshgrid and function values
	gg,pp = np.meshgrid(np.linspace(0,1.5, 100), np.linspace(1.5,0,100))
	z_array = [g_0(gg,pp), g_1(gg,pp), g_2(gg,pp)]

	#create figure 
	fig, axes = plt.subplots(1,3, constrained_layout=True)
	fig.suptitle('Coupling Strength Geometry', fontsize = 25)
	for n,ax in enumerate(axes,0):
		colormap = ax.imshow(z_array[n], extent=[0,1.5,0,1.5])
		ax.set_title(r'$|\Delta m| = ' + str(n)+'$')
		ax.set_xlabel(r'$\gamma$')
		ax.set_ylabel(r'$\phi$')
	fig.colorbar(colormap, ax=axes, location='bottom')

	plt.show()