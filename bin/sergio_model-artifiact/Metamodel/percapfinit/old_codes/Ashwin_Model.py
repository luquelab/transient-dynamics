
#Import libraries
#++++++++++++++++++++++++++++++++++++++++++++                                                       
import seaborn as sns
import numpy as np                                                                                  
from scipy.integrate import odeint                                                                  
import matplotlib.gridspec as gridspec                                                              
import matplotlib                                                                                   
import matplotlib.pyplot as plt                                                                     
from decimal import Decimal                                                                         
import math                                                                                         
from scipy.integrate import solve_ivp                                                               
import pandas as pd                                                                                 
from matplotlib.pyplot import figure                                                                
import copy                                                                                         
from Metamodel_Functions import *                               
import sys
#++++++++++++++++++++++++++++++++++++++++++++ 

#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++
#0.0 Weather model equation
#======================================
def Weather_model(t,Temp):

    return [a*(-Temp**4 + b_mu*Temp**2 - d_mu)]
#======================================

#++++++++++++++++++++++++++++++++++++++++++++ 

#Parameter values
I0=1366 #W/m2
sigma=5.6704e-8 #W/m2 K4
c=1e8 #kg K/s
e_sa=0.62
mu=1
#mu=0.0
b_2=1.69e-5 #1/K2
a_2=1.6927

#Derived parameters
a=e_sa*sigma/c
b_mu=mu*I0*b_2/(4*e_sa*sigma)
d_mu=-(mu*I0*(1-a_2))/(4*e_sa*sigma)

Temp_0=[290] #K
time_0=0   #years
time_f=300 #years
tau=100    #years

#1.3. Build time vector
#==============================
Resolution=100;Steps=time_f*Resolution + 1                                                          
time=np.linspace(time_0,time_f,Steps)                                                               
step=0.01                                          
#==============================

sol_model=solve_ivp(Weather_model,[time[0],time[-1]],Temp_0,method='RK45',dense_output=True,max_step=step)

z = sol_model.sol(time)

print(z)

plt.plot(time,z.T,color='gray',linewidth=1,label='Full')

plt.ylim([270,310])

plt.show()
                    

