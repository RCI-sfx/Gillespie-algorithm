#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:08:34 2020

@author: izzymillett
"""
from numpy.random import rand
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
k0 = 0.2 # transcription rate
k1 = 0.01 # degradation rate
Omega= 1. # cell size

mu = k0/k1

stoch = [+1, -1]
def propensities(x):
    return np.array((Omega*k0, k1*x))
def reactionTimes(x):
    a = propensities(x)
    aInv = [1/s if s>0 else np.inf for s in a]
    return -np.log(rand(2))*aInv

def ssa(x,tout):
    t=0.
    while t<tout:
        rt = reactionTimes(x)
        tau = np.min(rt)
        idx = np.argmin(rt)
        if t+tau>tout:
            t=tout
        else:
            t+=tau
            x+=stoch[idx]
    return x

z=1
k=0

x=k
dt=z
a=[]
for i in range(1000):
    x = ssa(x,dt)
    a.append(x)    
mean= np.mean(a[300:])
var= np.var(a[300:])
fano= (var/mean)
print(mean)
print(var)
print(fano)


x=k
dt=z
b=[]
for i in range(1000):
    x = ssa(x,dt)
    b.append(x) 

x=k
dt=z
c=[]
for i in range(1000):
    x = ssa(x,dt)
    c.append(x)
    
x=k
dt=z
d=[]
for i in range(1000):
    x = ssa(x,dt)
    d.append(x)

x=k
dt=z
f=[]
for i in range(1000):
    x = ssa(x,dt)
    f.append(x)

mean= np.mean(a[300:])
var= np.var(a[300:])
fano= (var/mean)
print(mean)
print(var)
print(fano)


mean= np.mean(b[300:])
var= np.var(b[300:])
fano= (var/mean)
print(mean)
print(var)
print(fano)

mean= np.mean(c[300:])
var= np.var(c[300:])
fano= (var/mean)
print(mean)
print(var)
print(fano)

mean= np.mean(d[300:])
var= np.var(d[300:])
fano= (var/mean)
print(mean)
print(var)
print(fano)

mean= np.mean(f[300:])
var= np.var(f[300:])
fano= (var/mean)
print(mean)
print(var)
print(fano)

plt.plot(a, label = "1 sim")
plt.plot(b, label = "2 sim")
plt.plot(c, label = "3 sim")
plt.plot(d, label = "4 sim")
plt.plot(f, label = "5 sim")

plt.xlabel('Time/seconds')
# Set the y axis label of the current axis.
plt.ylabel('[mRNA]')
# Set a title of the current axes.
plt.title('Gillespie simulation of [mRNA] ')
# show a legend on the plot
plt.legend()
# Display a figure.
plt.show()

