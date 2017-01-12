# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 14:51:19 2015

@author: hstrey
"""

from GaussianModels import g_n,vol1,vol2,modelFCS
import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
from scipy.interpolate import interp1d

# 3-d Gaussian focus FCS model
def g2(t,C,wxy,D):
    wz=1/wxy/wxy/np.pi**1.5/C
    return 1.0/(1+4*D*t/wxy/wxy)/np.sqrt(1+4*D*t/wz/wz)

modelFCS2 = lm.Model(g2,independent_vars=['t'])

def g2msd(msd,C,wxy):
    wz=1/wxy/wxy/np.pi**1.5/C
    return 1.0/(1+msd/wxy/wxy)/np.sqrt(1+msd/wz/wz)

# parameters that are assumed to be correctly describing the focus
w0=0.30
r0=0.17
lambdaex=0.488
lambdaem=0.517
n=1.33
a=1.044
D=420.0

# calculate g(t)-1 from these parameters
time=np.logspace(-7,-1,1000)
logtime=np.linspace(-7,-1,1000)

msd=time*D
logmsd=np.log10(msd)

# experimental noise spectrum
pf= [ -4.42459969e-04,4.96315309e-03,1.28443577e-01,1.50178679e-01,-2.76724306e+00]
p=np.poly1d(pf)
fitNoise=10**p(logtime)
real_g=g_n(time,D,1.0,w0,a,r0,lambdaex,lambdaem,n)

v1=vol1(a,r0,lambdaem,n)
v2=vol2(w0,a,r0,lambdaex,lambdaem,n)

print v1*v1/v2

# Concentration molecules per micrometer cubed
C=1/v1/v1*v2

real_g=(real_g-1.0)*v1*v1*6.022e-1/v2


resultG2 = modelFCS2.fit(real_g,t=time,
                   wxy=0.3,
                   C=lm.Parameter(value=C,vary =False),
                    D=lm.Parameter(value=D, vary = False),
                    weights=1.0/fitNoise)

print resultG2.fit_report()

plt.figure()
plt.semilogx(time,real_g)
plt.semilogx(time,resultG2.best_fit)

plt.figure()
plt.semilogx(time,resultG2.residual)

wxy=resultG2.values['wxy']
wz=1/wxy/wxy/np.pi**1.5/C

print wxy,wz

resultG = modelFCS.fit(real_g+1.0,t=time,
                       wxy=0.3,
                       wz=3.0,
                       C=lm.Parameter(value=C/6.022e-1,vary =False),
                       D=lm.Parameter(value=D, vary = False),
                       weights=1.0/fitNoise)

print resultG.fit_report()

plt.figure()
plt.semilogx(time,real_g+1.0)
plt.semilogx(time,resultG.best_fit)

plt.figure()
plt.semilogx(time,resultG.residual)

#now lets see how to get the mean square displacement out of the measurement
ginter=np.array([g2msd(d2,C,wxy) for d2 in msd])
# do the logarithmic interpolation
msdinter=interp1d(ginter,logmsd,bounds_error=False)

msdLog=msdinter(real_g)
msd_real=10**msdLog

#remove all nan
logtime=logtime[~np.isnan(msdLog)]
msdLog=msdLog[~np.isnan(msdLog)]

p=np.polyfit(logtime,msdLog,5)
p2=np.polyder(p)
msdDer=np.poly1d(p2)(logtime)

plt.figure()
plt.plot(logtime,msdDer)

