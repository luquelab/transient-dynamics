#26/5. A model for Sara's experiment. Question: given a single phage and a logistic growth curve for a bacterial culture, what's the concentration at which you would get the maximum number of phages?
#The model has three differential equations: one for sensitive bacteria, another one for infected bacteria, and another for phages. The equation for the infected bacteria provides a lag for phages so that burst is not instantaneous

#Import libraries
#++++++++++++++++++++++++++++++++++++++++++++  
import numpy as np   
from scipy.integrate import odeint   
import matplotlib.gridspec as gridspec
import matplotlib    
import matplotlib.pyplot as plt 
import seaborn as sns  
from decimal import Decimal 
import math            
from scipy.integrate import solve_ivp 
import pandas as pd 
from matplotlib.pyplot import figure
import copy  
from Metamodel_Functions import *
from matplotlib.colors import ListedColormap
import matplotlib.ticker as mticker
#++++++++++++++++++++++++++++++++++++++++++++

#0. FUNCTIONS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++
def Lotka_Volterra_Logistic(t,y):
    return [r*(1-y[0]/K)*y[0] - d*y[0]*y[1],
           Eclipse_Time*c*d*y[0]*y[1] - m*y[1]]

#Breaking condition for sensitive bottom 
#:::::::::::::::::::::::::::::::::::::::::::::
def Min_Volume_B(t,y): 
    return y[0] - 1  
Min_Volume_B.terminal=True  
Min_Volume_B.direction = 1 
#:::::::::::::::::::::::::::::::::::::::::::::

def Lagged_Lotka_Volterra_Logistic(t,y):
    return [r*(1-(y[0]+y[1])/K)*y[0] - d*(1-(y[0]+y[1])/K)*y[0]*y[2],
            d*(1-(y[0]+y[1])/K)*y[0]*y[2] - mu_p*y[1],
            c*mu_p*y[1] - m*y[2]]

#Avatar colormap
gyr=['#0F2347','#1C3F6E','#2E67A0','#5AACCF','#EFFC93','#80C271']
sns.palplot(sns.color_palette(gyr))
Avatar_cmap=ListedColormap(sns.color_palette(gyr))
#+++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++
#Parameters
Doubling_Time=1/3 #h
r=math.log(2)/Doubling_Time #h-1
K=1e8;d=3e-8
c=200
mu_p=1/(0.8*Doubling_Time)

Halving_Time_Ph=12 #h
m=math.log(2)/Halving_Time_Ph #h-1

#time
time_0=0; time_f=10
Resolution=100;Steps=time_f*Resolution + 1
time=np.linspace(time_0,time_f,Steps)
step=0.01

#Model solving
Solutions_L_V=[]
#Min concentrations
B0=[2,5,10,1e2,1e3,1e4];
#Max concentrations
#B0=[5e6,1e7,5e7,7e7,8e7,9e7]
P0=1

for concentration in B0:
    Events_Model=[Min_Volume_B]
    y0=[concentration,0,P0]
    sol_LV=solve_ivp(Lagged_Lotka_Volterra_Logistic,[time[0],time[-1]],y0,method='RK45',dense_output=True,events=Events_Model,max_step=step)

    z=sol_LV.sol(time)
    Solutions_L_V.append(z)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++

#Figure 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
cmap='Dark2';cmap_pieces=matplotlib.cm.get_cmap(cmap)
colors=\
[cmap_pieces(0.0),cmap_pieces(0.15),cmap_pieces(0.3),cmap_pieces(0.45),cmap_pieces(0.6),cmap_pieces(0.75), cmap_pieces(0.9), cmap_pieces(1.0)]
#Size of figures
cm = 1/2.54  # centimeters in inches
Width=8*cm;Height=5.0*cm #Width and height of plots 
#Fontsizes
size_axis=7;size_ticks=6;size_title=5


#Plots
fig=figure(figsize=(Width, Height), dpi=300)
width_line=1

cnt=0
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
fmt = mticker.FuncFormatter(g)

for solution in Solutions_L_V:
    plt.plot(time, solution[2].T,linewidth=width_line,color=colors[cnt],label="B0 = {}".format(fmt(B0[cnt])))
    cnt+=1

plt.axhline(y=K, color='k', linewidth=width_line,linestyle='dashed')
plt.axhline(y=1, color='k', linewidth=width_line,linestyle='dashed')
scale='linear'
if scale=='log':
    plt.yscale('log')
    plt.ylim(0.1, 1e12)
    
plt.xlabel('Time (h)',fontsize=size_axis)
plt.ylabel('Concentration (1/ml)',fontsize=size_axis)
plt.title('Phage concentration for initial bacterial concentrations',fontsize=size_axis)

plt.xticks(fontsize=size_ticks);
xticks_pos=[0, 0.1, 0.2,0.3, 0.4 ,0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
plt.yticks(xticks_pos,fontsize=size_ticks)
plt.ylim(0,1.1)
plt.legend(loc='best',fontsize=size_ticks,frameon=False) 
plt.show()
