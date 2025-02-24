#2021_10_04. Do a single experiment

#Import libraries
#+++++++++++++++++++++++++++++++++++++++
import numpy as np
from scipy.integrate import odeint
import math
from scipy.integrate import solve_ivp
import pandas as pd
from decimal import Decimal
from Metamodel_Functions import *
import pickle
from scipy import stats
import json
#++++++++++++++++++++++++++++++++

np.random.seed(1111)

#0. Functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Experiment dictionary
#==========================================================
def Build_Experiment_Dict(Number_of_Experiments, Measure_Interval_var):

    Measuring_Sample=5
    Measurements={}
    for Measure in range(Measuring_Sample):
        Measurements[Measure]=list(np.random.choice(\
                    Measure_Interval_var,4,replace=False))

    Experiments_dict_var={}
    for Experiment_Number in range(Number_of_Experiments):

        Measurements_test={}
        for Measure in range(Measuring_Sample):
            Measurements_test[Measure]=list(np.random.choice(Measure_Interval_var,4,replace=False))

        Experiments_dict_var[Experiment_Number]=Measurements_test

    return Experiments_dict_var
#==========================================================

#0.2. Statistical Tests
#============================================================================== 
def Statistical_Tests(Sample1, Sample2):
    
    Variance_1=Sample1.var(ddof=0)
    Variance_2=Sample2.var(ddof=0)
    
    if Variance_1!=Variance_2:
        Ttest=stats.ttest_ind(Sample1, Sample2, equal_var=False)

    else:
        Ttest=stats.ttest_ind(Sample1, Sample2)
        
    WMWtest=stats.mannwhitneyu(Sample1, Sample2)
    
    return Ttest, WMWtest
#==============================================================================

#Extract VMR dataframe
#==============================================================================
def Extract_VMR_fun(Snapshot_df_var):
    Cells_df=Snapshot_df_var[Snapshot_df_var["Bioagent"]=='CFU']
    Cells_df.reset_index(drop=True, inplace=True)
    
    Phages_df=Snapshot_df_var[Snapshot_df_var["Bioagent"]=='PFU']
    Phages_df.reset_index(drop=True, inplace=True) 

    VMR=Phages_df['Concentration']/Cells_df['Concentration']
    VMR_series = pd.Series(VMR, name="Concentration")
    VMR_df_var=VMR_series.to_frame()
    VMR_df_var['Experiment']=Cells_df['Experiment']
    VMR_df_var['Measurement']=Cells_df['Measurement']
    VMR_df_var['Time']=Cells_df['Time'];VMR_df_var['Bioagent']='VMR'
    VMR_df_var['Type']=Cells_df['Type']
    VMR_df_var['Total_Type']=VMR_df_var['Bioagent'] + VMR_df_var['Type']

    return VMR_df_var
#==============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. Main
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Input_Path='/home/sergio/work/Github/needle-finder/data/Data_Marisa/'
Output_Path='/home/sergio/work/Github/needle-finder/results/Data/Results_Marisa/'

Medium='LB'
Data=Input_Path + 'Data_Marisa_' + Medium + '.csv'
Raw_Data=pd.read_csv(Data,sep='\t',index_col='Time (hr)') 

#Initial Conditions
#===================================
B0_WT=Raw_Data['WT'][0.0]
B0_sp=Raw_Data['spoT-'][0.0]

B0=B0_WT;T0=0
Initial_Conditions=[B0,T0]
#===================================

#Time and resolution
#======================================
time_0=0;time_f=35
Resolution=100
Steps=time_f*Resolution+1
time=np.linspace(time_0,time_f,Steps)
step=0.01 
#====================================== 

#Parameters
#========================================================================
r_wt=0.471; r_sp=0.382
K_wt=Raw_Data['WT'][30.0]; K_sp=Raw_Data['spoT-'][30.0]
mu_wt=3.9935398580806127e-11; mu_sp=8.481314602001232e-12
c=125
m=0.003

Parameters_wt={'r':r_wt, 'K':K_wt, 'mu':mu_wt, 'c':c, 'm': m}
Parameters_sp={'r':r_sp, 'K':K_sp, 'mu':mu_sp, 'c':c, 'm': m}
#========================================================================

#Model Solution
#========================================================================
Solution_wt=\
Solve_Experiment_Induction(Parameters_wt, Initial_Conditions, time,step)

Solution_sp=\
Solve_Experiment_Induction(Parameters_sp, Initial_Conditions, time,step)

#Save solutions to dictionary and then to dataframe
#---------------------------------------------------------------------
Solutions_Dict={\
    'wt':{'CFU':Solution_wt.sol(time)[0],'PFU':Solution_wt.sol(time)[1]},
    'sp':{'CFU':Solution_sp.sol(time)[0],'PFU':Solution_sp.sol(time)[1]}}

Solutions_df=Solutions_to_Dataframe(Solutions_Dict,time)
print(Solutions_df)
#---------------------------------------------------------------------

#========================================================================

#Measure Time: theoretical and experimental
#============================================================================
#Randomly choose where to measure from the pool
Measurement_Time=10
Measure_Err=0.5 #h

#The threshold is an error of 50%
PFU_Error=10; PFU_Threshold=PFU_Error*2

#Extract dataframes with Experimental Data and pvalues
#=================================================================================
Measurement_Interval=Measurement_Interval_fun(Measurement_Time,Measure_Err,time)

Number_of_Experiments=1
Experiments=Build_Experiment_Dict(Number_of_Experiments,Measurement_Interval)

#Save experiment times and expected concentrations to dataframe                
#---------------------------------------------------------------------
Snapshot_Experiments_df=Snapshots_to_Dataframe(\
Experiments,Solutions_Dict,time,Measurement_Time)
#---------------------------------------------------------------------
#=================================================================================
        
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Save data
# Snapshot_Experiments_df.to_csv(Output_Path + 'Snapshot_Experiments.csv')

# Solutions_df.to_csv(Output_Path + 'Solutions.csv')

# Whole_Space_df.to_csv(Output_Path + 'Solutions' + '_' + str(Measurement_Times) + '_Measurements_' + str(Number_of_Experiments) + '_Number_of_Experiments' + '.csv')

# Pvalues_df.to_csv(Output_Path + 'Pvalues' + '_' + str(Measurement_Times) + '_Measurements_' + str(Number_of_Experiments) + '_Number_of_Experiments' + '.csv')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

