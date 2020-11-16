import numpy as np
import matplotlib.pyplot as plt

def dicke_regime_flops_red(t, rabi, eta, n_av):
    r = 0
    for n in range(500):
        r += ((n_av)/(n_av + 1))**(n+1)*(np.sin(rabi*eta*t*np.sqrt(n + 1)/2)**2)/(n_av + 1)
    return r

def dicke_regime_flops_blue(t, rabi, eta, n_av):
    r = 0
    for n in range(500):
        r += ((n_av)/(n_av + 1))**(n)*(np.sin(rabi*eta*t*np.sqrt(n + 1)/2)**2)/(n_av + 1)
    return r

def dicke_regime_flops_carrier(t, rabi, eta, n_av):
    r = 0
    for n in range(500):
        r += ((n_av)/(n_av + 1))**(n)*(np.sin(rabi*t*(1 - n*eta**2)/2)**2)/(n_av + 1)
    return r

if __name__ == "__main__":

    t = np.linspace(0,1000,1000)
    #in microseconds  

    #Now let out pi time be 30us (i.e. how well you can identify the carrier)
    Rabi = np.pi/50
    #Rabi now in MHz 

    #corresponds to axial frequency of 1MHz
    eta = 0.068

    n_av = 0.5

    plt.plot(t, dicke_regime_flops_blue(t, Rabi, eta, n_av))
    plt.plot(t, dicke_regime_flops_red(t, Rabi, eta, n_av))
    plt.plot(t, dicke_regime_flops_carrier(t, Rabi, eta, n_av))
    plt.show()