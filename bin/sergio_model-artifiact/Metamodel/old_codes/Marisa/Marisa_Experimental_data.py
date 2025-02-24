#19/08/2021. Script to analyse Marisa's data in CFUs. The script does 3 things:
#1. Reads Marisas data and converts to dataframe and calculate empirical rates. You have to choose between rich media (LB) and minimal media (MM)
#2. Solve induction model numerically
#3. Plot results

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


#1. Read data and infer growth rates (exponential and linear) & induction rates
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Input_Path='/home/sergio/work/Github/needle-finder/data/Data_Marisa/'
Medium='LB' #LB or MM
Raw_Data=Input_Path + 'Data_Marisa_' + Medium + '.csv'
Data=pd.read_csv(Raw_Data,sep='\t',index_col='Time (hr)')
print(Data)

#Calculate exponential and linear growth rates
#=================================================================
if Medium=='LB':
    t0=0;t_CX=5.5;t_f=30
    K=1.4e9
elif Medium=='MM':
    t0=0;t_CX=4;t_f=24
    K=1.4e8

Exponential_Growth_Rates=Growth_Rates(t_CX,t0, Data)
print(Exponential_Growth_Rates)

Linear_Growth_Rates=Growth_Rates(t_f,t_CX, Data)
print(Linear_Growth_Rates)
#=================================================================

#Calculate induction Rates
#=================================================================
c=125
if Medium=='LB':
    Phages={'spoT-': 4e1,'WT': 1.7e2,'spoT- +CX':8e3, 'WT +CX':2.2e4}
elif Medium=='MM':
    Phages={'spoT-':2.1e3 ,'WT':5.6e3 ,'spoT- +CX':5e3 , 'WT +CX':2.5e4 }

    
Terminal_Concentrations=Data.iloc[-1]
print("Terminal concentrations")
print(Terminal_Concentrations)
Induction_Rates={}
for key in Phages:
    Induced_Bacteria=Phages[key]/c #phages divided by burst size
    Induced_Ratio=Induced_Bacteria/(Induced_Bacteria+Terminal_Concentrations[key])
    Induction_Rates[key]=Induced_Ratio/(t_f-t_CX)
    print(key)
    print(Phages[key])
    print(Terminal_Concentrations[key])
    print(Induced_Bacteria)
    print(Induced_Ratio)
    print(t_CX)

print("Induction rates")
print(Induction_Rates)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#Model Simulation
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.2.1. Time vector for exponential phase
#:::::::::::::::::::::::::::::::::::::::::::::::
step=0.01
Simulation_t0=t0;Simulation_t_CX=t_CX
t_exp=np.linspace(Simulation_t0,Simulation_t_CX,101)
#:::::::::::::::::::::::::::::::::::::::::::::::

#1.2.2. Solve model in Exponential phase
#:::::::::::::::::::::::::::::::::::::::::::::::
Solutions_Exp={}
Initial_Concentrations=Data.loc[t0]

for key in Phages:
    r=Exponential_Growth_Rates[key]
    print(key)
    print(Exponential_Growth_Rates[key])
    mu_i=Induction_Rates[key]
    
    L0_exp=Initial_Concentrations[key];T0_exp=1
    y0_exp=[L0_exp,T0_exp]

    solution=solve_ivp(Induction_Model,[t_exp[0],t_exp[-1]],y0_exp,method='RK45',dense_output=True,events=None,max_step=step)
    #Save solution in dictionary
    Solutions_Exp[key]=solution
#:::::::::::::::::::::::::::::::::::::::::::::::

#1.2.2. Time vector for linear phase
#:::::::::::::::::::::::::::::::::::::::::::::::
step=0.01;Simulation_tf=t_f
t=np.linspace(Simulation_t_CX,Simulation_tf,101)
#:::::::::::::::::::::::::::::::::::::::::::::::

#1.2.3. Solve model - linear phase
#:::::::::::::::::::::::::::::::::::::::::::::::
Solutions={}

CX_Initial_Concentrations=Data.loc[Simulation_t_CX]

for key in Phages:
    
    r=Linear_Growth_Rates[key]
    mu_i=Induction_Rates[key]
    L0=CX_Initial_Concentrations[key]
    T0=1

    y0_1=[L0,T0]
    solution=solve_ivp(Induction_Model,[t[0],t[-1]],y0_1,method='RK45',dense_output=True,events=None,max_step=step)
    
    Solutions[key]=solution
#:::::::::::::::::::::::::::::::::::::::::::::::
#=================================================================================

#Print data
print(Data.index.tolist())

for key in Solutions:
    
    Lexp=Solutions_Exp[key].y[0]
    print(Lexp[-1])
    print(type(Lexp))
    L=Solutions[key].y[0]
    print(L[0])
    Total_L=np.concatenate((Lexp,L), axis=None)
    print(len(Lexp))
    print(len(L))
    print(len(Total_L))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
#2. PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fig=plt.figure(figsize=(13, 4))
gs=gridspec.GridSpec(2, 1)
gs.update(left=0.06,right=0.96,bottom=0.12,top=0.92,wspace=0.2,hspace=0.45)

colors=['orange','b']
styles=['solid','dashed']
plot_style={'WT':['orange','solid','o'], 'WT +CX':['orange', 'dashed','^'],\
            'spoT-':['b', 'solid','o'],  'spoT- +CX': ['b', 'dashed','^']}

ax00=plt.subplot(gs[0,0])
for key in Solutions:

    Lexp=Solutions_Exp[key].y[0]
    t_exp=np.linspace(Simulation_t0,Simulation_t_CX,len(Lexp))
    plt.plot(t_exp,Lexp,linestyle=plot_style[key][1], color=plot_style[key][0])
    
    L=Solutions[key].y[0]

    t_plot=np.linspace(Simulation_t_CX,Simulation_tf,len(L))

    plt.plot(t_plot,L,linestyle=plot_style[key][1],color=plot_style[key][0],label=key)
    plt.scatter(Data.index.tolist(), Data[key],color=plot_style[key][0],marker=plot_style[key][2],label=key)

    
plt.ylabel('Concentration (cells/ml)')
plt.xlabel('time (h)')
plt.title('CFUs')

    
plt.legend(loc='upper left')


ax00=plt.subplot(gs[1,0])

for key in Solutions:

    T=Solutions[key].y[1]

    plt.plot(t_plot,T,linestyle=plot_style[key][1], color=plot_style[key][0],label=key)

#    plt.yscale('log')

plt.ylabel('Concentration (Phages/ml)')
plt.xlabel('time (h)')
plt.title('Phages')
plt.legend(loc='upper left')
plt.show()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
