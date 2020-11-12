import numpy as np 
import matplotlib.pyplot as plt 

f1 = 1
f2 = 1.2

#fibre noise model 
phi_fibre = lambda t: 0.01*t**2

t = np.linspace(0,200,10000)
x = lambda t: np.cos(f1*t)
y = lambda t: np.cos(f2*t + phi_fibre(t))



plt.plot(t,x(t)+y(t))
plt.show() 
