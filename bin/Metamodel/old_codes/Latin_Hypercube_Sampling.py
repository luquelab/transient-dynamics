
#Import libraries
#+++++++++++++++++++++++++++++++++++++++++++++++++++
from Metamodel_Function import *
import numpy as np
from random import shuffle
import matplotlib.gridspec as gridspec
import time
import matplotlib.pyplot as plt
import seaborn as sns
#+++++++++++++++++++++++++++++++++++++++++++++++++++


#1. Parameter values and intervals
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Sampling Points
N=10

#Growth rates
range_r=[0.00295, 0.00704]

#Infection rates             
range_d=[7.2e-10, 3.7e-7]

#Burst size
range_c=[50,600]

#Phage decay rate
range_m=[0.04,0.2363]

#Carrying capacity
range_K=[1e6,1e8]

Ranges_Parameters={'r':range_r,'d':range_d,'c':range_c,'m':range_m,'K':range_K}

Samples_Parameters=Latin_Hypercube_Sampling(Ranges_Parameters, N, 3333)
print(Samples_Parameters)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3. Lotka-Volterra Solution
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Initial Conditions
#================
B0=1000
T0=10000
Initial_Conditions=[B0,T0]
#================

#Time and resolution
#======================================
time_0=0
time_f=10
Resolution=1
Steps=time_f*Resolution+1
time=np.linspace(time_0,time_f,Steps)
step=1
#======================================

counter_coexistence=0
counter_bacteria=0

B_eq_vector=[]
T_eq_vector=[]
VMR_eq_vector=[]

K_vector=[];r_vector=[]
for sampling in Samples_Parameters:
    
    Parameters=Samples_Parameters[sampling]

    r=Parameters['r'];K=Parameters['K'];m=Parameters['m'];
    c=Parameters['c'];d=Parameters['d']

    B_eq=m/(c*d)
    T_eq=(r/d)*(1-m/(c*d*K))
    VMR_eq=T_eq/B_eq
    
    if K<(m/(c*d)):
        counter_bacteria+=1
    elif K>(m/(c*d)):
        counter_coexistence+=1
        B_eq_vector.append(B_eq)
        T_eq_vector.append(T_eq)
        VMR_eq_vector.append(VMR_eq)


    
    K_vector.append(K)
    r_vector.append(r)    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print(counter_bacteria)
print(counter_coexistence)

#4. Plots
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fig=plt.figure(figsize=(15, 5))
gs=gridspec.GridSpec(1, 3)
gs.update(left=0.07,right=0.95,bottom=0.08,top=0.94,wspace=0.2,hspace=0.2)

ax_00=plt.subplot(gs[0,0])
sns.histplot(B_eq_vector,log_scale=True)
plt.title('Bacterial Concentration')

ax_01=plt.subplot(gs[0,1])
sns.histplot(T_eq_vector, log_scale=True)
plt.title('Phage Concentration')

ax_02=plt.subplot(gs[0,2])
sns.histplot(VMR_eq_vector, log_scale=True)
plt.title('VMR')

plt.show()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
