# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 11:32:43 2015

@author: hstrey
"""

import numpy as np
from scipy.integrate import nquad,quad
import matplotlib.pyplot as plt
from matplotlib import rcParams
from GaussianModels import vol1,vol2,vol2c

rcParams.update({'figure.autolayout': True})

w1=0.2
w2=0.3
r1=0.123
r2=0.1
lambdaex1=0.488
lambdaem1=0.519
lambdaex2=0.633
lambdaem2=0.657
n=1.33
a=0.522

def RR(z):
    return r0*r0+(lambdaem*z/np.pi/r0/n)**2
    
def k(z):
    return 1-np.exp(-2*a*a/(r0*r0+(lambdaem*z/np.pi/r0/n)**2))

def cefz_int(r,z):
    return 2*np.greater(RR(z),r**2)*r/RR(z)

def cef_int(r,alpha,r_prime,R):
    return np.greater(R**2,r**2+r_prime**2-2*r*r_prime*np.cos(alpha))*r
    
def cef_z_geo(z):
    return np.where(RR(z)>a**2,a**2/RR(z),1.0)
    
z=np.linspace(0,5,50)
cef=[quad(cefz_int,0,a,args=(zz)) for zz in z]
cefg=cef_z_geo(z)
kk=k(z)

#r_prime=np.linspace(0,1.0,50)
#cef_r=[nquad(cef_int,[[0,0.5],[0,2*np.pi]],args=(r_p,0.170)) for r_p in r_prime]

plt.figure(figsize=(4,3))
#plt.plot(z, cef)
plt.plot(z,kk)
plt.plot(z,cefg)
plt.xlabel("z in micrometers")
plt.ylabel("Transmission")
#plt.figure()
#plt.plot(r_prime,cef_r)
