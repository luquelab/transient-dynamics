#19/08/2021. Script to analyse Marisa's data in CFUs
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
from matplotlib.pyplot import figure
from Metamodel_Functions import *
from lmfit import minimize, Parameters, Parameter, report_fit

#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.0. Calculate Growth rates
#=================================================================================
def Growth_Rates_fun(final_t, initial_t, Dataframe):
    
    Final_Concentration=Dataframe.loc[final_t]
    Initial_Concentration=Dataframe.loc[initial_t]
    
    Log_Diff_Concentration=\
    np.log(Final_Concentration) - np.log(Initial_Concentration)
    
    Time=final_t - initial_t
    Growth_Rates=Log_Diff_Concentration/Time
    
    return Growth_Rates
#=================================================================================

#0.1. Calculate Induction Rates
#=================================================================================
def Induction_Rates_fun(PFUs_end,Terminal_Concentrations_var,Burst_Size,t_0,t_f):
    Induction_Rates_var={}

    for key in Phages:
        Induced_Bacteria=Phages[key]/Burst_Size
        Induced_Ratio=\
        Induced_Bacteria/(Induced_Bacteria+Terminal_Concentrations_var[key])
        Induction_Rates_var[key]=Induced_Ratio/(t_f-t_CX)

    return Induction_Rates_var
#=================================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Input_Path='/home/sergio/work/Github/needle-finder/data/Data_Marisa/'
Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/'

Medium='MM'
Raw_Data=Input_Path + 'Data_Marisa_' + Medium + '.csv'
Data=pd.read_csv(Raw_Data,sep='\t',index_col='Time (hr)')

#Growth rates, carrying capacities, and Induction rates
#=================================================================
c=125;m=0.003

if Medium=='LB':
    t0=0;t_CX=5.5;t_f=30

    #Initial concentrations
    B0_wt=Data['WT'][0.0]
    B0_sp=Data['spoT-'][0.0]
    T0=0
    #Carrying capacity
    K_wt=Data['WT'][t_f]; K_sp=Data['spoT-'][t_f]
    #Growth rates
    Growth_Rates=Growth_Rates_fun(t_CX,t0,Data)
    r_wt=Growth_Rates['WT']
    r_wt_cx=Growth_Rates['WT +CX']
    r_sp=Growth_Rates['spoT-']
    r_sp_cx=Growth_Rates['spoT- +CX']
    #Possibly right PFUs at the end of the experiment
#    Phages={'spoT-': 9e3,'WT': 6e4,'spoT- +CX':1e4, 'WT +CX':9e4}
    #possibly wrong PFUs
    Phages={'spoT-': 4e1,'WT': 1.7e2,'spoT- +CX':8e3, 'WT +CX':2.2e4}
    Terminal_Concentrations=Data.iloc[-1]
    Induction_Rates=Induction_Rates_fun(Phages,Terminal_Concentrations,c,t_CX,t_f)
    mu_wt=Induction_Rates['WT'];mu_wt_cx=Induction_Rates['WT +CX']
    mu_sp=Induction_Rates['spoT-'];mu_sp_cx=Induction_Rates['spoT- +CX']
    
elif Medium=='MM':
    t0=0;t_CX=4;t_f=24

    #Initial concentrations
    B0_wt=Data['WT'][0.0]
    B0_sp=Data['spoT-'][0.0]
    T0=0
    #Carrying capacity
    K_wt=Data['WT'][t_f]; K_sp=Data['spoT-'][t_f]
    #Growth rates
    #Growth rates of LB that make the model much more accurate
    r_wt_LB=0.4707509373925001
    r_sp_LB=0.38182924161501314    
    Growth_Rates=Growth_Rates_fun(t_CX,t0,Data)
    r_wt=Growth_Rates['WT'];r_wt_cx=Growth_Rates['WT +CX']
    r_sp=Growth_Rates['spoT-'];r_sp_cx=Growth_Rates['spoT- +CX']
    #Induction rates
    Phages={'spoT-': 4e1,'WT': 1.7e2,'spoT- +CX':8e3, 'WT +CX':2.2e4}
    Terminal_Concentrations=Data.iloc[-1]
    Induction_Rates=Induction_Rates_fun(Phages,Terminal_Concentrations,c,t_CX,t_f)
    mu_wt=Induction_Rates['WT'];mu_wt_cx=Induction_Rates['WT +CX']
    mu_sp=Induction_Rates['spoT-'];mu_sp_cx=Induction_Rates['spoT- +CX']
#=================================================================

#1.1. Model
#=================================================================

#Initial Conditions
#--------------------------------
Initial_Conditions_wt=[B0_wt,T0]
Initial_Conditions_sp=[B0_sp,T0]
#--------------------------------

#Time and resolution
#---------------------------------------
Resolution=100
Steps=t_f*Resolution+1
time=np.linspace(t0,t_f,Steps)
step=0.001
#---------------------------------------

#Parameters
#--------------------------------------------------------------
Parameters_wt={'r':r_wt, 'K':K_wt, 'mu':mu_wt, 'c':c, 'm': m}
Parameters_sp={'r':r_sp, 'K':K_sp, 'mu':mu_sp, 'c':c, 'm': m}
#--------------------------------------------------------------

#Model Solution
#----------------------------------------------------------------------------
Solution_wt=Solve_Experiment_Induction(\
Parameters_wt, Initial_Conditions_wt, time,step)           

Solution_sp=Solve_Experiment_Induction(\
Parameters_sp, Initial_Conditions_sp, time,step)

#Save solutions to dictionary and then to dataframe
#.........................................................................
Solutions_Dict={\
    'wt':{'CFU':Solution_wt.sol(time)[0],'PFU':Solution_wt.sol(time)[1]},
    'sp':{'CFU':Solution_sp.sol(time)[0],'PFU':Solution_sp.sol(time)[1]}}
Solutions_df=Solutions_to_Dataframe(Solutions_Dict,time)
Solutions_df['Total_Type']=Solutions_df['Bioagent'] + Solutions_df['Type']
#.........................................................................

#----------------------------------------------------------------------------

#=================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
colors=['orange','b']
styles=['solid','dashed']
plot_style={'WT':['orange','solid','o'], 'WT +CX':['orange', 'dashed','^'],\
            'spoT-':['b', 'solid','o'],  'spoT- +CX': ['b', 'dashed','^']}


#Parameters and fontsizes
#========================
size_axis=16
size_ticks=14
size_title=18
#======================== 

figure(figsize=(17, 6), dpi=80)
#Plot CFUs
#=================================================================================
for column in Data:

    plt.scatter(Data.index.tolist(), Data[column], color=plot_style[column][0],marker=plot_style[column][2],label=column, s=100)

    plt.plot(Data.index.tolist(), Data[column], color=plot_style[column][0], linestyle=plot_style[column][1],linewidth=3)

Solutions_CFU_wt=Solutions_df[Solutions_df["Total_Type"]=="CFUwt"]
Solutions_CFU_sp=Solutions_df[Solutions_df["Total_Type"]=="CFUsp"]

plt.plot(Solutions_CFU_wt["Time"], Solutions_CFU_wt["Concentration"], color=plot_style['WT'][0], linestyle='dotted', linewidth=3, label='WT model')

plt.plot(Solutions_CFU_sp["Time"], Solutions_CFU_sp["Concentration"], color=plot_style['spoT-'][0], linestyle='dotted', linewidth=3, label='spoT- model')
    
plt.ylabel('Concentration (cells/ml)',fontsize=size_axis)
plt.yscale('log')
plt.xlabel('time (h)',fontsize=size_axis)
#plt.title('CFUs in %s Media'  % Medium, fontsize=size_title)
plt.legend(loc='upper left')
plt.xticks(fontsize=size_ticks)
plt.yticks(fontsize=size_ticks)

plt.ylim(ymax = 2.5e9, ymin = 6e6)

plt.savefig(Output_Path+'Experiment_CFUs_in_'+str(Medium)+'_Media.png',dpi=300)
plt.show()
#=================================================================================

#Plot PFUs
#=================================================================================
figure(figsize=(17, 6), dpi=80)
Solutions_PFU_wt=Solutions_df[Solutions_df["Total_Type"]=="PFUwt"]
Solutions_PFU_sp=Solutions_df[Solutions_df["Total_Type"]=="PFUsp"]

plt.plot(Solutions_PFU_wt["Time"], Solutions_PFU_wt["Concentration"], color=plot_style['WT'][0], linestyle='dotted', linewidth=3,label='wild type')
plt.plot(Solutions_PFU_sp["Time"], Solutions_PFU_sp["Concentration"], color=plot_style['spoT-'][0], linestyle='dotted', linewidth=3,label='spoT-')
    
plt.ylabel('Concentration (cells/ml)',fontsize=size_axis)
plt.yscale('log')
plt.xlabel('time (h)',fontsize=size_axis)
#plt.title('PFUs in %s Media'  % Medium, fontsize=size_title)
plt.legend(loc='upper left')

plt.xticks(fontsize=size_ticks)
plt.yticks(fontsize=size_ticks)


plt.ylim(ymax = 5000, ymin = 0)
plt.axhline(y=1)

plt.savefig(Output_Path+'Experiment_PFUs_in_'+str(Medium)+'_Media.png',dpi=300)

plt.show()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
