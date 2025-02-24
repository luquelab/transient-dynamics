import numpy as np
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd
from Metamodel_Functions import *

from matplotlib.pyplot import figure


#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    
    for key in PFUs_end:
        Induced_Bacteria=PFUs_end[key]/Burst_Size
        Induced_Ratio=\
        Induced_Bacteria/(Induced_Bacteria+Terminal_Concentrations_var[key])
        Induction_Rates_var[key]=Induced_Ratio/(t_f-t_0)
        
    return Induction_Rates_var
#=================================================================================

#0.2. Extract data from cvs file
#=================================================================================
def Extract_Parameters_fun(Data_var, Medium_var):    
    c=125;m=0.003    
    if Medium_var=='LB':
        t0=0;t_CX=5.5;t_f=30
        #Initial concentrations
        B0_wt=Data['WT'][0.0];B0_sp=Data['spoT-'][0.0];T0=1e-8
        #Carrying capacity
        K_wt=Data['WT'][t_f];K_sp=Data['spoT-'][t_f]
        #Growth rates
        Growth_Rates=Growth_Rates_fun(t_CX,t0,Data)
        r_wt=Growth_Rates['WT'];r_wt_cx=Growth_Rates['WT +CX']
        r_sp=Growth_Rates['spoT-'];r_sp_cx=Growth_Rates['spoT- +CX']
        #Induction rates
        Phages={'spoT-': 4e1,'WT': 1.7e2,'spoT- +CX':8e3, 'WT +CX':2.2e4}
        Terminal_Concentrations=Data.iloc[-1]
        Induction_Rates=Induction_Rates_fun(\
        Phages,Terminal_Concentrations,c,t_CX,t_f)
        mu_wt=Induction_Rates['WT'];mu_wt_cx=Induction_Rates['WT +CX']
        mu_sp=Induction_Rates['spoT-'];mu_sp_cx=Induction_Rates['spoT- +CX']

    elif Medium_var=='MM':
        t0=0;t_CX=4;t_f=24
        #Initial concentrations
        B0_wt=Data['WT'][0.0];B0_sp=Data['spoT-'][0.0];T0=1e-8
        #Carrying capacity
        K_wt=Data['WT'][t_f]; K_sp=Data['spoT-'][t_f]
        #Growth rates
        Growth_Rates=Growth_Rates_fun(t_CX,t0,Data)
        r_wt=Growth_Rates['WT'];r_wt_cx=Growth_Rates['WT +CX']
        r_sp=Growth_Rates['spoT-'];r_sp_cx=Growth_Rates['spoT- +CX']
        #Induction rates
#        Phages={'spoT-': 4e1,'WT': 1.7e2,'spoT- +CX':8e3, 'WT +CX':2.2e4}
        Phages={'spoT-': 9e3,'WT': 6e4,'spoT- +CX':1e4, 'WT +CX':9e4}
        Terminal_Concentrations=Data.iloc[-1]
        Induction_Rates=\
        Induction_Rates_fun(Phages,Terminal_Concentrations,c,t_CX,t_f)
        mu_wt=Induction_Rates['WT'];mu_wt_cx=Induction_Rates['WT +CX']
        mu_sp=Induction_Rates['spoT-'];mu_sp_cx=Induction_Rates['spoT- +CX']


    Params_wt={   'r':r_wt,'K':K_wt,'mu':mu_wt,   'c':c,'m':m,'B0':B0_wt,'T0':T0}
    Params_wt_cx={'r':r_wt,'K':K_wt,'mu':mu_wt_cx,'c':c,'m':m,'B0':B0_wt,'T0':T0}
    Params_sp={   'r':r_sp,'K':K_sp,'mu':mu_sp,   'c':c,'m':m,'B0':B0_sp,'T0':T0}
    Params_sp_cx={'r':r_sp,'K':K_sp,'mu':mu_sp_cx,'c':c,'m':m,'B0':B0_sp,'T0':T0}

    Params={'WT':Params_wt, 'WT +CX':Params_wt_cx, 'spoT-':Params_sp, 'spoT- +CX':Params_sp_cx}
    time_points={'t0':t0,'t_CX':t_CX,'t_f':t_f}

    return Params, time_points
#=================================================================================

#===========================================
def residual(ps, ts, data):
    Initial_Conditions = ps['B0'].value, ps['T0'].value
    model=Solve_Experiment_Induction(ps, Initial_Conditions, ts, step)
    model=model.sol(ts)
    model=model.T
    return (model - data).ravel()
#===========================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Induction Model
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Time and initial conditions
#=================================
# Length=100
# time = np.linspace(0, 100, Length)
#=================================

#Parameter values
#======================================================================
Input_Path='/home/sergio/work/Github/needle-finder/data/Data_Marisa/'
Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/'
Medium='MM'
Strain='WT'
Raw_Data=Input_Path + 'Data_Marisa_' + Medium + '.csv'
Data=pd.read_csv(Raw_Data,sep='\t',index_col='Time (hr)')
#======================================================================

#Growth rates, carrying capacities, and Induction rates
#=================================================================

Parameters_Dict,time_points=Extract_Parameters_fun(Data, Medium)

Parameters_Strain=Parameters_Dict[Strain]
Experimental_Data_raw=Data[Strain]

t0=time_points['t0']
t_f=time_points['t_f']
#=================================================================

#Time and initial conditions
#===============================
Length=t_f
time = np.linspace(t0, t_f, Length)
Initial_Conditions=[ Parameters_Strain['B0'], Parameters_Strain['T0'] ]
#===============================

#Solve model
#===============================
step=0.1
model_wt_raw=Solve_Experiment_Induction(\
Parameters_Strain,Initial_Conditions, time,step)
model_wt=model_wt_raw.sol(time)
model_wt=model_wt.T
#===============================

#Add noise
Noise_wt=np.random.normal(size=model_wt.shape)
data_wt=model_wt

#Choose data from the model at random
#=====================================================================
Time_indices=np.linspace(0,Length-1,Length)
Random_Indices=np.random.choice(Time_indices, 10, replace=False)
Random_Indices=[int(random_index) for random_index in Random_Indices]
Random_Indices=np.sort(Random_Indices)
Small_data=data_wt[Random_Indices]
Small_time=time[Random_Indices]
#=====================================================================

#Experimental data and pseudo experimental data (Phages)
#=====================================================================
Experimental_CFU=Experimental_Data_raw.to_list()
Experimental_times=Experimental_Data_raw.index.to_list()
Pseudo_exp_PFU=[]
for Experimental_t in Experimental_times:
    time_index=min(range(len(time)), key=lambda i: abs(time[i]-Experimental_t))
    Pseudo_exp_PFU.append(model_wt[time_index][1])

Experimental_Data=[]
for element in zip(Experimental_CFU,Pseudo_exp_PFU):
    Experimental_Data.append(element)

Experimental_Data=np.asarray(Experimental_Data)
#=====================================================================

#Set parameters
#===========================================================
print(Parameters_Strain['r'])
for Parameter in Parameters_Strain:
    r=Parameters_Strain['r']
    K=Parameters_Strain['K']
    mu=Parameters_Strain['mu']
    c=Parameters_Strain['c']
    m=Parameters_Strain['m']
    B0=Parameters_Strain['B0']
    T0=Parameters_Strain['T0']
    

params_induction = Parameters()
params_induction.add('B0', value= float(B0), min=0.5*B0, max=1.5*B0)
params_induction.add('T0', value=float(T0), min=0.5*T0, max=1.5*T0)
params_induction.add('r', value=r, min=0.5*r, max=1.5*r)
params_induction.add('K', value=K, min=0.5*K, max=1.5*K)
params_induction.add('mu', value=mu, min=0.5*mu, max=1.5*mu)
params_induction.add('c', value=c, min=0.5*c, max=1.5*c)
params_induction.add('m', value=m, min=0.5*m, max=1.5*m)
#===========================================================

# fit model and find predicted values
#=====================================================================
result=minimize(residual,params_induction,args=(Experimental_times, Experimental_Data),method='leastsq')
final = Experimental_Data + result.residual.reshape(Experimental_Data.shape)
#=====================================================================

#Extract values for plot
#========================================
Residuals=result.residual
Residuals_CFU=Residuals[::2]
Residuals_PFU=Residuals[1:][::2]

CFUs_Model=model_wt[:,0]
PFUs_Model=model_wt[:,1]

CFUs_Exp=Experimental_Data[:,0]
PFUs_Exp=Experimental_Data[:,1]

CFUs_Fit=final[:,0]
PFUs_Fit=final[:,1]

print(result.params)
#========================================


#PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Parameters and fontsizes
#========================
size_axis=14
size_ticks=12
size_title=16
#========================  

#================================================================================
figure(figsize=(17, 6), dpi=80)
plt.plot(Experimental_times,Residuals_CFU,'o',color='blue',linewidth=3,label='cells')
plt.axhline(y=0, color='k')

#Axes, title and ticks
plt.ylabel('Residuals',fontsize=size_axis)
plt.xlabel('time (h)',fontsize=size_axis)
plt.title('Residuals of '+str(Strain) + ' in ' + str(Medium), fontsize=size_title)
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)

#axes limits
plt.ylim(-1.5e7,1.3e7)

plt.savefig(Output_Path+'Model_Fit_Residuals_CFUs_'+str(Strain)+'_' + str(Medium)+ '.png',dpi=300)

plt.show()

#================================================================================

#================================================================================
figure(figsize=(17, 6), dpi=80)
plt.plot(Experimental_times,Residuals_PFU,'o',color='red',linewidth=3,label='phages')                                         
plt.axhline(y=0, color='k')

#Axes, title and ticks
plt.ylabel('Residuals',fontsize=size_axis)
plt.xlabel('time (h)',fontsize=size_axis)
plt.title('Residuals of '+str(Strain) + ' in ' + str(Medium), fontsize=size_title)
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)

#Legend and scale                                               
#plt.yscale('log')
#plt.legend(loc='upper left')
plt.savefig(Output_Path+'Model_Fit_Residuals_PFUs_'+str(Strain)+'_' + str(Medium)+ '.png',dpi=300)

plt.show()


#================================================================================

#================================================================================
figure(figsize=(17, 6), dpi=80)

plt.plot(time, CFUs_Model, linestyle='dotted', color='blue', linewidth=3,label='Original model') 
plt.plot(Experimental_times, CFUs_Exp, 'o', color='blue', label='experimental data')
plt.plot(Experimental_times, CFUs_Fit, '-', color='blue', linewidth=3, label='Fitted model');

plt.ylabel('Concentration (cells/ml)',fontsize=size_axis)
plt.yscale('log')

plt.xlabel('time (h)',fontsize=size_axis) 
plt.title('Model vs. Data for ' + str(Strain) + ' in ' + str(Medium) + ' media',fontsize=size_title)
plt.legend(loc='lower right')
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)

#axes limits
plt.ylim(6.7e6, 1.5e8)

plt.savefig(Output_Path+'Model_Fit_CFUs_'+str(Strain)+'_' + str(Medium) + '.png',dpi=300)

plt.show()
#================================================================================

#================================================================================
figure(figsize=(17, 6), dpi=80)
plt.plot(time, PFUs_Model, linestyle='dotted', color='red',linewidth=3,label='model')
plt.plot(Experimental_times, PFUs_Exp, 'o', color='red', label='experimental_data')
plt.plot(Experimental_times, PFUs_Fit, '-', color='red',  linewidth=2, label='model fitting')
plt.axhline(y=0, color='k')

plt.ylabel('Concentration (phages/ml)',fontsize=size_axis)
plt.yscale('log')
plt.xlabel('time (h)',fontsize=size_axis)
plt.title('Model vs. Data for ' + str(Strain) + ' in ' + str(Medium) + ' media',fontsize=size_title)
plt.legend(loc='lower right')
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)
plt.savefig(Output_Path+'Model_Fit_PFUs_'+str(Strain)+'_' + str(Medium)+ '.png',dpi=300)
plt.show() 
# display fitted statistics
report_fit(result)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
