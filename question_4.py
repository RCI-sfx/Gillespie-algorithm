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
div_index = []
for i in range(10000):
   k0=random.uniform(0.2,0.4) # this can make K0 randomonly distrubted
   r1=k0*Omega # rate for mRNA production
   r2=k1*m #rate for mRNA degrdation
   rtot= r1 +r2 #total reaction rates 
   randtot=random.rand()*rtot #  random number between 0-1 times by rtot
   if randtot <= r1: #determining which reaction to fire 
       m = m+1
   else: #Rxn 2 fires (mRNA degradation)
       m = m-1
   randtime=random.exponential(1/rtot) #working out time interval   
   t= t + randtime 
    
   if   t >= 1200*divtime: #working out when to divide
        divtime = divtime + 1 
        m = np.random.binomial(m, 0.5) #allow odd division of integers +/- 0.5 equally 
        div_index.append(i) #store information 
    
   m_store.append(m)
   t_store.append(t)
  
x= t_store[div_index[0]:div_index[1]] #working out one cell division mRNA 

mean= np.mean(x)
var= np.var(x)
fano= (var/mean)
print(mean)
print(var)
print(fano)

mean= np.mean(m_store) #working out the total one cell division mRNA 
var= np.var(m_store)
fano= (var/mean)
print(mean)
print(var)
print(fano)

plt.plot(t_store, m_store, color='orchid')
plt.xlabel('Time/seconds')
# Set the y axis label of the current axis.
plt.ylabel('[mRNA]')
# Set a title of the current axes.
plt.title('Gillespie simulation of [mRNA] in cell division ')
# show a legend on the plot
# Display a figure.
plt.show()