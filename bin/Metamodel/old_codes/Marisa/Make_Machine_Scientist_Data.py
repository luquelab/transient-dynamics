#27/10/2021. This script will solve a simple L-V system and save data in a datrame
#for the machine scientist


#Import libraries
#+++++++++++++++++++++++++++++++++++
from Metamodel_Functions import *
import pandas as pd
import numpy as np
#+++++++++++++++++++++++++++++++++++


#1. Main
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Output_Path='/home/sergio/work/Github/needle-finder/results/Data/Results_Marisa/'

#========================================================================
#Fix parameters and initial conditions
r=0.9;K=2e6;d=3e-8;c=150;m=2.8e-3
Parameters={'r':r,'K':K,'d':d,'c':c,'m':m}

B0=1e3;T0=1e4
Initial_Conditions=[B0,T0]

#Time vector
t0=0;tf=35;
Resolution=100
Steps=tf*Resolution+1
time=np.linspace(t0,tf,Steps)
step=0.1

#Solve
Solution=Solve_Lotka_Volterra(Parameters,Initial_Conditions,time,step)
#========================================================================

#Put solutions and derivatives into a dataframe
#========================================================================
Bacteria=Solution.y[0]
Phage=Solution.y[1]

Bacteria_Derivative_1=[(Bacteria[i1]-Bacteria[i1-1])/step for i1 in range(1,len(Bacteria))]
Bacteria_Derivative_2=\
[r*Bacteria[i2]-d*Phage[i2]*Bacteria[i2] for i2 in range(1,len(Bacteria))]


Phage_Derivative_1=[(Phage[j1]-Phage[j1-1])/step for j1 in range(1,len(Phage))]
Phage_Derivative_2=\
[c*d*Phage[j2]*Bacteria[j2] - m*Phage[j2] for j2 in range(1, len(Phage))]

#Put data in dict
Solutions_Dict={'B':Bacteria[1:], 'T':Phage[1:], \
                'dB/dt_1':Bacteria_Derivative_1,'dT/dt_1':Phage_Derivative_1,
                'dB/dt_2':Bacteria_Derivative_2,'dT/dt_2':Phage_Derivative_2}

#Put dict into df
Data_machine_scientist_df=pd.DataFrame.from_dict(Solutions_Dict)
#========================================================================



#Save data
#========================================================================
Data_machine_scientist_df.to_csv(Output_Path + 'Data_machine_scientist.csv')
#========================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
