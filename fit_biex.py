import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math 
import os
import csv 


somedata1=np.loadtxt(fname="av_rot_corr.dat")
#somedata2=np.loadtxt(fname="SET.2/stay_connect_time_corr")
#somedata3=np.loadtxt(fname="SET.3/stay_connect_time_corr")
#y1=somedata1[:,1]
#y2=somedata2[:,1]
#y3=somedata3[:,1]
x=somedata1[:,0]
#y=(y1+y2+y3)/3
y=somedata1[:,1]
#print(y)
# Function to calculate the exponential with constants a and b
def exponential(x, a0, a1, a2):
    return (a0*np.exp(-a1*x))+ ((1-a0)*np.exp(-(a2*x)))

# Fit the dummy exponential data
pars, cov = curve_fit(f=exponential, xdata=x, ydata=y, p0=[0.2, 0.1,0.01], bounds=(0, 1))
print(pars)
print(cov)
fig=plt.figure()
# Plot the fit data as an overlay on the scatter data
plt.scatter(x, y, s=20, color='#00b3b3', label='Data')
plt.plot(x, exponential(x, *pars), linestyle='--', linewidth=2, color='black')
#ax.set_xlim(0.0, 200)
#ax.set_ylim(0.0, 1)
plt.show()

# Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
stdevs = np.sqrt(np.diag(cov))# Calculate the residuals
res = y - exponential(x, *pars)

ts=(pars[0]*(1/pars[1])) + ((1-pars[0])*(1/pars[2]))
print(type(ts))

f = open('ts.dat', 'w')
f.write("%12.6f\n" % (ts))
f.close()
