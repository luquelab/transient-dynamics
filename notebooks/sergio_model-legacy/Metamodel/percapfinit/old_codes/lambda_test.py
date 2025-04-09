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
import matplotlib.ticker as mtick                                                                   
#++++++++++++++++++++++++++++++++++++++++++++


def Weather_model_kick(t,lambda_var):
    return [rho*lambda_var*(1-lambda_var)]



rho=0.19
time_0=0   #years                                                                                 
time_f=150 #years 

y0=[10**(-6)]

#1.3. Build time vector    
#==============================                                                                     
Resolution=100;Steps=time_f*Resolution + 1                                                          
time=np.linspace(time_0,time_f,Steps)                                                               
step=0.01                                                                                           
#==============================


sol_lambda=solve_ivp(Weather_model_kick,[time[0],time[-1]],y0,method='RK45',dense_output=True,max_step=step)

z=sol_lambda.sol(time)

lambda_function=z.T
b_2_0=1.69e-5;b_2_f=1.835e-5
b2=(1-lambda_function)*b_2_0 + lambda_function*b_2_f

size_axis=7;size_ticks=6;size_title=5
cm = 1/2.54  # centimeters in inches
Width=8*cm;Height=4*cm #Width and height of plots

plt.plot(time,lambda_function,color='k',linewidth=2,label='Equilibrium model')
#plt.plot(time,b2,color='gray',linewidth=2,label='Equilibrium model')
plt.xlabel('Time (years)',fontsize=size_axis)
plt.ylabel('$\lambda$ ',fontsize=size_axis)
#plt.ylabel('b_2 (K$\mathrm{^{-2}}$)',fontsize=size_axis)
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)


plt.show()
