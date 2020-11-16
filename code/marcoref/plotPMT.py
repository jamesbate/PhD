data_dir_trics = 'C:/Users/James/OneDrive - OnTheHub - The University of Oxford/phd/data/'
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 16})

def read_data(file):
    data = np.genfromtxt(file, delimiter="\t",usecols=range(0,106),skip_header=1)
    return data[:,6:], data[:,2]

def calculate_probability(n,t1,t2,tot):
	p0 = 0
	p1 = 0
	p2 = 0
	for i in range(len(n)):
		if i<t1:
			p0+=n[i]
		elif i>=t1 and i<t2:
			p1+=n[i]
		else:
			p2+=n[i]
	return p0/tot,p1/tot,p2/tot


def poissonian(x,mean,A):
	return A*(mean**x*np.exp(-mean))/(factorial(x))

def cos(x,A,phi):
	frequency = 1 #N/2 
	return A*np.cos(2*np.pi*frequency*x + phi)


def histograms(x,length,treshold1,treshold2):
	plt.figure()

	x1 = np.array([])
	for j in range(length):
		x1 = np.concatenate((x1,x[j,:]))

	n, bins, patches = plt.hist(x1, np.arange(-0.5,100.5,1), facecolor='g', alpha=0.75)

	plt.axvline(x=treshold1, color='r')
	plt.axvline(x=treshold2, color='r')

	bin_center = [ (bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]

	popt, pcov = curve_fit(poissonian, bin_center[:treshold1], n[:treshold1])

	print(popt)

	xdata = np.linspace(-1,treshold1,1000)
	plt.plot(xdata, poissonian(xdata, *popt),'b--')


	popt, pcov = curve_fit(poissonian, bin_center[treshold1:treshold2], n[treshold1:treshold2],p0=[(treshold1+treshold2/2),10])

	print(popt)

	xdata = np.linspace(treshold1,treshold2,1000)
	plt.plot(xdata, poissonian(xdata, *popt),'b--')


	popt, pcov = curve_fit(poissonian, bin_center[treshold2:], n[treshold2:],p0=[2*treshold2,10])

	print(popt)

	xdata = np.linspace(treshold2,100,1000)
	plt.plot(xdata, poissonian(xdata, *popt),'b--')

	plt.grid(True)

folder = "181333" #pop 181542 #par 181333

x,phase = read_data(data_dir_trics+folder+"/PMT1_2.txt")
length,cycles = np.shape(x)

treshold1 = 15
treshold2 = 45

histograms(x,length,treshold1,treshold2)

p0 = np.zeros(length)
p1 = np.zeros(length)
p2 = np.zeros(length)

dp0 = np.zeros(length)
dp1 = np.zeros(length)
dp2 = np.zeros(length)

for i in range(length):
	n,bins = np.histogram(x[i,:], np.arange(-0.5,100.5,1))
	p0[i],p1[i],p2[i] = calculate_probability(n,treshold1,treshold2,cycles)
	dp0[i] = np.max([np.sqrt(p0[i]*(1-p0[i])/cycles),1/(cycles+2)])
	dp1[i] = np.max([np.sqrt(p1[i]*(1-p1[i])/cycles),1/(cycles+2)])
	dp2[i] = np.max([np.sqrt(p2[i]*(1-p2[i])/cycles),1/(cycles+2)])


plt.figure()
plt.errorbar(phase,p0,yerr = dp0,fmt="--o",label="p0")
plt.errorbar(phase,p1,yerr = dp1,fmt="--o",label="p1")
plt.errorbar(phase,p2,yerr = dp2,fmt="--o",label="p2")

plt.ylabel("Probability")
plt.xlabel(r"Phase/$2\pi$")
plt.grid(True)
plt.legend()


plt.figure()


parity = p0+p2-p1
dparity = np.sqrt(dp0**2+dp1**2+dp2**2)

plt.errorbar(phase,parity,yerr = dparity,fmt="--o")
plt.xlabel(r"Phase/$2\pi$")
plt.ylabel(r"Parity $p_0+p_2-p_1$")
plt.grid(True)


popt, pcov = curve_fit(cos, phase,parity,sigma= dparity, p0=[0.9,0],absolute_sigma=False) # Absolute sigma??
perr = np.sqrt(np.diag(pcov))
print("parity upper")
print(popt,perr)

a1 = popt[0]
da1 = perr[0]

xdata = np.linspace(0,2,1000)
y = cos(xdata, *popt)

upper = y[y>0]
xupper = xdata[y>0]
#plt.plot(xdata, cos(xdata, *popt),'r--')

#coherence = np.abs(popt[0])
#dcoherence = perr[0]

#stupidate

#plt.figure()
parity2 = parity[parity<0]
phase2 = phase[parity<0]
dparity2 = dparity[parity<0]

#plt.errorbar(phase2,parity2,yerr = dparity2,fmt="--o")
#plt.xlabel(r"Phase/$2\pi$")
#plt.ylabel(r"Parity $p_0+p_2-p_1$")
#plt.grid(True)

popt, pcov = curve_fit(cos, phase2,parity2,sigma= dparity2, p0=[0.9,0],absolute_sigma=False) # Absolute sigma??
perr = np.sqrt(np.diag(pcov))
print("parity lower")
print(popt,perr)
a2 = popt[0]
da2 = perr[0]

y =  cos(xdata, *popt)
lower = y[y<0]
xlower = xdata[y<0]

plt.plot(xupper,upper,'r--')
plt.plot(xlower,lower,'r--')

coherence = np.abs(a1+a2)/2
dcoherence = np.sqrt(da1**2 + da2**2)/2


#################################################### POPULATION

folder = "181542" #pop 181542 #par 181333

x,phase = read_data(data_dir_trics+folder+"/PMT1_2.txt")
length,cycles = np.shape(x)


treshold2 = 41

histograms(x,length,treshold1,treshold2)


p0 = np.zeros(length)
p1 = np.zeros(length)
p2 = np.zeros(length)

dp0 = np.zeros(length)
dp1 = np.zeros(length)
dp2 = np.zeros(length)

for i in range(length):
	n,bins = np.histogram(x[i,:], np.arange(-0.5,100.5,1))
	p0[i],p1[i],p2[i] = calculate_probability(n,treshold1,treshold2,cycles)
	dp0[i] = np.max([np.sqrt(p0[i]*(1-p0[i])/cycles),1/(cycles+2)])
	dp1[i] = np.max([np.sqrt(p1[i]*(1-p1[i])/cycles),1/(cycles+2)])
	dp2[i] = np.max([np.sqrt(p2[i]*(1-p2[i])/cycles),1/(cycles+2)])


plt.figure()
plt.errorbar(phase,p0,yerr = dp0,fmt="--o",label="p0")
plt.errorbar(phase,p1,yerr = dp1,fmt="--o",label="p1")
plt.errorbar(phase,p2,yerr = dp2,fmt="--o",label="p2")

plt.ylabel("Probability")
plt.xlabel(r"Random stuff")
plt.grid(True)
plt.legend()


p0mean = np.average(p0,weights = 1/dp0**2)
dp0mean = np.sqrt(1/np.sum(1/dp0**2))
p2mean = np.average(p2,weights = 1/dp2**2)
dp2mean = np.sqrt(1/np.sum(1/dp2**2))

p1mean = np.average(p1,weights = 1/dp1**2)
dp1mean = np.sqrt(1/np.sum(1/dp1**2))

print(p0mean,dp0mean)

print(p2mean,dp2mean)

print(p1mean,dp1mean)


#plt.figure()
#plt.plot(phase,p1,"-o")
#plt.axhline(y=p1mean)

plt.figure()

plt.hist(p1, np.arange(-0.005,0.105,0.01), facecolor='g', alpha=0.75)
plt.axvline(x=np.mean(p1))

print("std:",p1mean,np.std(p1),np.std(p1)/np.sqrt(len(p1)))

############# FIDELITY


fidelity = (p0mean + p2mean + coherence)/2
dfidelity = 0.5*np.sqrt(dp0mean**2+dp2mean**2 +dcoherence**2)
print("Coherence",coherence,dcoherence)
print("Fidelity",fidelity,dfidelity)


plt.show()

