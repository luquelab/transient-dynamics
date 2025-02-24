#21/09. First model for Marisa's project

#Import libraries
#+++++++++++++++++++++++++++++++++++++++
import numpy as np
from scipy.integrate import odeint
import math
from scipy.integrate import solve_ivp
import pandas as pd
from decimal import Decimal
from Metamodel_Functions import *
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
    
    Experiments_dict_var={}
    
    for Experiment_Number in range(Number_of_Experiments):
        Measurements_test={}
        for Measure in range(Measuring_Sample):
            Measurements_test[Measure]=list(np.random.choice(Measure_Interval_var,4,replace=False))

        Experiments_dict_var[Experiment_Number]=Measurements_test

    return Experiments_dict_var
#==========================================================

#0.2. Solutions to dataframe
#==========================================================
def Solutions_to_Dataframe(Solutions_Dict_var):
    
    List_df=[]
    for key1 in Solutions_Dict_var:
        for key2 in Solutions_Dict_var[key1]:
            for snapshot in zip(Solutions_Dict_var[key1][key2],time):
                List_df.append([snapshot[1],snapshot[0],key2,key1])

    Solutions_df_var=pd.DataFrame(List_df,columns=[\
                "Time","Concentration", "Bioagent","Type"])

    return Solutions_df_var
#==========================================================


#0.3. Solutions LHS dataframe
#==========================================================
def Solutions_LHS_to_Dataframe(Solutions_LHS_Dict_var):
    
    List_df=[]
    for key0 in Solutions_LHS_Dict_var:
        for key1 in Solutions_LHS_Dict_var[key0]:
            for key2 in Solutions_LHS_Dict_var[key0][key1]:
                for snapshot in zip(Solutions_LHS_Dict_var[key0][key1][key2],time):
                    List_df.append([snapshot[1],snapshot[0],key2,key1,key0])

    Solutions_LHS_df_var=pd.DataFrame(List_df,columns=[\
    "Time","Concentration", "Bioagent","Type","Sampling"])

    return Solutions_LHS_df_var
#==========================================================

#0.3. Experiments to Dataframe
#================================================================================
def Snapshots_to_Dataframe(Experiments_Dict_var, Solutions_Dict_var, time_var, Measure_time_var):

    List_df=[]

    for Experiment_Number in Experiments_Dict_var:
        for Measurement in Experiments_Dict_var[Experiment_Number]:
            
            #Save snapshots of the model to list for Dataframe
            Times_Snapshot=Experiments_Dict_var[Experiment_Number][Measurement]

            counter_time=0
            for key1 in Solutions_Dict_var:
                for key2 in Solutions_Dict_var[key1]:

                    time_experiment=Times_Snapshot[counter_time]

                    index=np.where(time_var==time_experiment)

                    List_df.append([Measure_time_var,Experiment_Number,Measurement,time_var[index][0],Solutions_Dict_var[key1][key2][index][0],key2,key1])

                    counter_time+=1

    Snapshot_df_var=pd.DataFrame(List_df,columns=\
    ["Theoretic Time","Experiment","Measurement","Time","Concentration","Bioagent","Type"])
    

    return Snapshot_df_var
#================================================================================

#0.4. Statistical Tests
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

#Latin Hypercube Sampling
#------------------------------------------------------------------
Interval_Width=0.1 #Ten percent of the numerical value
N=5#Sampling size

range_r_wt=[(1-Interval_Width)*r_wt, (1+Interval_Width)*r_wt]
range_r_sp=[(1-Interval_Width)*r_sp, (1+Interval_Width)*r_sp]

range_K_wt=[(1-Interval_Width)*K_wt, (1+Interval_Width)*K_wt]
range_K_sp=[(1-Interval_Width)*K_sp, (1+Interval_Width)*K_sp]

range_mu_wt=[(1-Interval_Width)*mu_wt, (1+Interval_Width)*mu_wt]
range_mu_sp=[(1-Interval_Width)*mu_sp, (1+Interval_Width)*mu_sp]

range_c=[(1-Interval_Width)*c, (1+Interval_Width)*c]
range_m=[(1-Interval_Width)*m, (1+Interval_Width)*m]

Ranges_wt_LHS={\
'r':range_r_wt,'K':range_K_wt,'mu':range_mu_wt,'c':range_c,'m': range_m}

Ranges_sp_LHS={\
'r':range_r_sp,'K':range_K_sp,'mu':range_mu_sp,'c':range_c,'m': range_m}


LHS_wt=LHS(Ranges_wt_LHS, 5, 1111)
LHS_sp=LHS(Ranges_sp_LHS, 5, 1111)
#------------------------------------------------------------------

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

Solutions_Model_df=Solutions_to_Dataframe(Solutions_Dict)
#---------------------------------------------------------------------


#Solution within LHS
#---------------------------------------------------------------------
Solutions_LHS={}
for Sampling in LHS_wt:
    
    Solution_wt_lhs=\
    Solve_Experiment_Induction(LHS_wt[Sampling], Initial_Conditions, time,step)

    Solution_sp_lhs=\
    Solve_Experiment_Induction(LHS_sp[Sampling], Initial_Conditions, time,step)
    
    Solutions_Dict_lhs={\
    'wt':{'CFU':Solution_wt_lhs.sol(time)[0],'PFU':Solution_wt_lhs.sol(time)[1]},
    'sp':{'CFU':Solution_sp_lhs.sol(time)[0],'PFU':Solution_sp_lhs.sol(time)[1]}}

    Solutions_LHS[Sampling]=Solutions_Dict_lhs

Solutions_Model_LHS_df=Solutions_LHS_to_Dataframe(Solutions_LHS)
#---------------------------------------------------------------------

#========================================================================

#Measure Time: theoretical and experimental
#============================================================================
#Choose the pool of times where measurements are possible
Time_Sampling_Space=np.linspace(time_0,time_f,time_f+1)

#Randomly choose where to measure from the pool
Sample_Size=36
Measurement_Times=np.random.choice(Time_Sampling_Space,Sample_Size,replace=False)
Measure_Err=0.5 #h

#The threshold is an error of 50%
PFU_Error=10; PFU_Threshold=PFU_Error*2

#Extract dataframes with Experimental Data and pvalues
#=================================================================================
Pvalues_list=[];Non_Excluded_Times=[]
for Measurement_Time in Measurement_Times:
    
    Measurement_Interval=Measurement_Interval_fun(\
    Measurement_Time,Measure_Err,time)

    #Choose number of Experiments per time of measurement
    #---------------------------------------------------------------------
    Number_of_Experiments=100
    Experiments=Build_Experiment_Dict(Number_of_Experiments,Measurement_Interval)
    #---------------------------------------------------------------------

    #Save experiment times and expected concentrations to dataframe
    #---------------------------------------------------------------------
    Snapshot_Experiments_df=Snapshots_to_Dataframe(\
    Experiments,Solutions_Dict,time,Measurement_Time)
    #---------------------------------------------------------------------

    #First condition of measurement
    #---------------------------------------------------------------------
    PFUs=Snapshot_Experiments_df[Snapshot_Experiments_df["Bioagent"]=='PFU']
    Non_Observable_PFUs=PFUs[PFUs["Concentration"]< PFU_Threshold]
    
    if len(Non_Observable_PFUs) <1:
        
        Non_Excluded_Times.append(Measurement_Time)
        try:
            Whole_Space_dataframe=\
            Whole_space_df.append(Snapshot_Experiments_df)
        except NameError:
            Whole_Space_df=Snapshot_Experiments_df
    #---------------------------------------------------------------------

    #Merge two columns
    #---------------------------------------------------------------------
    Snapshot_Experiments_df["Total_Type"]=\
    Snapshot_Experiments_df["Bioagent"] + Snapshot_Experiments_df["Type"]
    #---------------------------------------------------------------------

    #Extract VMR
    #...............................................
    VMR_df=Extract_VMR_fun(Snapshot_Experiments_df)
    VMR_wt=VMR_df[VMR_df['Total_Type']=='VMRwt']
    VMR_sp=VMR_df[VMR_df['Total_Type']=='VMRsp']
    #...............................................

    #Filter Snapshots by bioagent
    #...............................................
    Snapshots_CFU_wt=Snapshot_Experiments_df[\
    Snapshot_Experiments_df['Total_Type']=='CFUwt']
    
    Snapshots_CFU_sp=Snapshot_Experiments_df[\
    Snapshot_Experiments_df['Total_Type']=='CFUsp']
    
    Snapshots_PFU_wt=Snapshot_Experiments_df[\
    Snapshot_Experiments_df['Total_Type']=='PFUwt']
    
    Snapshots_PFU_sp=Snapshot_Experiments_df[\
    Snapshot_Experiments_df['Total_Type']=='PFUsp'] 
    #...............................................

    for Experiment in range(Number_of_Experiments):
        #Statistical tests on CFU_sp vs CFU_wt
        #-------------------------------------------------------------------------
        Snapshots_CFU_wt_Experiment=\
        Snapshots_CFU_wt[Snapshots_CFU_wt['Experiment']==Experiment]

        Snapshots_CFU_sp_Experiment=\
        Snapshots_CFU_sp[Snapshots_CFU_sp['Experiment']==Experiment]

        Ttest_Cells,WMWtest_Cells=Statistical_Tests(Snapshots_CFU_wt_Experiment['Concentration'],Snapshots_CFU_sp_Experiment['Concentration'])
        #-------------------------------------------------------------------------

        #Statistical tests on PFU_sp vs PFU_wt
        #-------------------------------------------------------------------------
        Snapshots_PFU_wt_Experiment=\
        Snapshots_PFU_wt[Snapshots_PFU_wt['Experiment']==Experiment]

        Snapshots_PFU_sp_Experiment=\
        Snapshots_PFU_sp[Snapshots_PFU_sp['Experiment']==Experiment]

        Ttest_Phages,WMWtest_Phages=Statistical_Tests(Snapshots_PFU_wt_Experiment['Concentration'],Snapshots_PFU_sp_Experiment['Concentration'])
        #-------------------------------------------------------------------------

        #Statistical tests on VMR_sp vs VMR_wt
        #-------------------------------------------------------------------------
        VMR_wt_Experiment=VMR_wt[VMR_wt['Experiment']==Experiment]
        VMR_sp_Experiment=VMR_sp[VMR_sp['Experiment']==Experiment]
        Ttest_VMRs,WMWtest_VMRs=Statistical_Tests(\
        VMR_wt_Experiment['Concentration'],VMR_sp_Experiment['Concentration'])
        #-------------------------------------------------------------------------

        #List of pvalues
        #-------------------------------------------------------------------------
        Pvalues_list.append([Measurement_Time,"Cells",\
                             WMWtest_Cells[1],Experiment])
        
        Pvalues_list.append([Measurement_Time,"Phages",\
                             WMWtest_Phages[1],Experiment])
        
        Pvalues_list.append([Measurement_Time,"VMR",\
                             WMWtest_VMRs[1],Experiment])
        #-------------------------------------------------------------------------
        
#=================================================================================
print("Non excluded times")
print(Non_Excluded_Times)
        
Pvalues_df=pd.DataFrame(Pvalues_list,columns=[\
        "Time","Bioagent", "Pvalues","Experiment"])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Save data
#==========================================================================
Solutions_Model_df.to_csv(Output_Path + 'Solutions.csv')

Solutions_Model_LHS_df.to_csv(Output_Path + 'Solutions_LHS_'+\
        str(Interval_Width)+'_Interval_variation.csv')

Whole_Space_df.to_csv(Output_Path + 'Solutions' + '_' + str(Sample_Size) + '_Measurements_' + str(Number_of_Experiments) + '_Number_of_Experiments' + '.csv')

Pvalues_df.to_csv(Output_Path + 'Pvalues' + '_' + str(Sample_Size) + '_Measurements_' + str(Number_of_Experiments) + '_Number_of_Experiments' + '.csv')
#==========================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

