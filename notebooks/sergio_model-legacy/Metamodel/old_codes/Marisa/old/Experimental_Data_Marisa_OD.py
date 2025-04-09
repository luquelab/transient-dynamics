#12/8/2021. Script to analyse Marisa's data

import pandas as pd
import math
import numpy as np
import xlrd
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Calculate growth rates
#=================================================================================
def Growth_Rates(final_t, initial_t, Dataframe):
    Final_Concentration=Dataframe.loc[final_t]
    Initial_Concentration=Dataframe.loc[initial_t]

    Log_Diff_Concentration=np.log(Final_Concentration) - np.log(Initial_Concentration)
    Time=final_t - initial_t
    
    Growth_Rate=Log_Diff_Concentration/Time

    return Growth_Rate
#=================================================================================


#ODE's model
#================================
def Induction_Model(t,y):

    L,T=y
    return [r*y[0] - mu_i*y[0],
            c*mu_i*y[1]   ]
#================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#1. MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Raw_Data=pd.read_csv('Data_Marisa.csv',sep='\t',index_col='Time (hr)')

#Conversion of OD to cells/ml for E. Coli OD600 of 1.0 = 8 x 108 cells/ml
Raw_Data[['spoT-','WT','spoT- +CX', 'WT +CX']]=\
Raw_Data[['spoT-','WT','spoT- +CX', 'WT +CX']]*8e8

#Growth rates
#=================================================================================
#Instantaneous
Instant_Growth_Rate=np.log(Raw_Data[['spoT-','WT','spoT- +CX', 'WT +CX']]).diff()/0.5

#Exponential and Linear
t0=1
t1=6
t2=21.5

Exponential_Growth_Rates=Growth_Rates(t1, t0, Raw_Data)
Linear_Growth_Rates=Growth_Rates(t2, t1, Raw_Data)
#=================================================================================

#Induction rate
#=================================================================================
c=125
Phages={'spoT-': 4e1,'WT': 2e2,'spoT- +CX':8e3, 'WT +CX':2e4}

print(Raw_Data)
Terminal_Concentrations=Raw_Data.iloc[-1]
print(Terminal_Concentrations['spoT-'])

Induction_Rates={}
for key in Phages:
    
    Induced_Bacteria=Phages[key]/c
    Induced_Ratio=Induced_Bacteria/(Induced_Bacteria+Terminal_Concentrations[key])
    Induction_Rates[key]=Induced_Ratio/(t2-t1)
#=================================================================================

#Model Simulation
#======================================================================================

#1.2.1. Time and step
#:::::::::::::::::::::::::::::::::::::::::::::::
step=0.01
Simulation_t0=1
Simulation_tf=30
t=np.linspace(Simulation_t0,Simulation_tf,101)
#:::::::::::::::::::::::::::::::::::::::::::::::

#1.2.2. Solve model
#:::::::::::::::::::::::::::::::::::::::::::::::
Solutions={}

Initial_Concentrations=Raw_Data.loc[t1]
print(Initial_Concentrations)

for key in Phages:

    r=Linear_Growth_Rates[key]
    mu_i=Induction_Rates[key]
    L0=Initial_Concentrations[key]
    T0=1

    y0_1=[L0,T0]
    print(L0)
    
    solution=solve_ivp(Induction_Model,[t[0],t[-1]],y0_1,method='RK45',dense_output=True,events=None,max_step=step)

    Solutions[key]=solution
#:::::::::::::::::::::::::::::::::::::::::::::::

#======================================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fig=plt.figure(figsize=(13, 4))

colors=['orange','b']
styles=['solid','dashed']
plot_style={'WT':['orange','solid'], 'WT +CX':['orange', 'dashed'],\
            'spoT-':['b', 'solid'],  'spoT- +CX': ['b', 'dashed']}

for key in Solutions:
    
    L=Solutions[key].y[0]
    T=Solutions[key].y[1]
                                                                     
    t_plot=np.linspace(1,t[-1],len(L)) 


    plt.plot(t_plot,L,linestyle=plot_style[key][1], color=plot_style[key][0],label=key)
    plt.plot(t_plot,T,color='r',label='Phages')
#    plt.yscale('log')

    plt.ylabel('Concentration (cells/ml)')
    plt.xlabel('time (h)')

plt.legend(loc='upper left')
plt.show()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

