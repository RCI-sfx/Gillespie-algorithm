#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 15:11:17 2020

@author: izzymillett
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 12:36:43 2020

@author: izzymillett
"""
from numpy import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
k0 = 0.2 # transcription rate
k1 = 0.01 # degradation rate
Omega= 1. # cell size
k2= 5 # protein translation
m=0  # starting mRNA
p=0 # starting protein 
t=0 # starting time 
divtime= 1

m_store = [m]
t_store= [t]
p_store= [p]
div_index = []
for i in range(10000000):
   r1=k0*Omega # rate for mRNA production
   r2=k1*m #rate for mRNA degrdation
   r3=m*k2 #
   rtot= r1 +r2 + r3 #random number
   randtot=random.rand()*rtot #  Rxn 1 fires (mRNA production)
   if randtot <= r1:
       m = m+1
   elif randtot >=r1 and randtot <r2+r1  : #Rxn 2 fires (mRNA degradation)
       m = m-1
   else:
       p=p+1
   randtime=random.exponential(1/rtot)      
   t= t + randtime 
 
   if t >= 1200 * divtime:
        divtime = divtime + 1 
        m = np.random.binomial(m, 0.5)
        p = np.random.binomial(p, 0.5)
        div_index.append(i)
    
   m_store.append(m)
   t_store.append(t)
   p_store.append(p)
   

mean= np.mean(p_store) #working out the total one cell division mRNA 
var= np.var(p_store)
print(mean)
print(var)
fano= (var/mean)
print(fano)


plt.plot(t_store, p_store)
plt.xlabel('Time/seconds')
# Set the y axis label of the current axis.
plt.ylabel('[Protein]')
# Set a title of the current axes.
plt.title('Gillespie simulation of [Protien]')

x=p_store[div_index[0]:div_index[1]] #working out one cell division mRNA 
mean= np.mean(x)
var= np.var(x)
print(mean)
print(var)
fano= (var/mean)
print(fano)
     
del p_store[0:div_index[1]]
mean= np.mean(m_store) #working out the total one cell division mRNA 
var= np.var(m_store)
fano= (var/mean)
print(mean)
print(var)
print(fano)
