#23/05/2022. Considering m=r dominant timescales and starting
#from equilibrium concentrations, I run the model and then I
#kick viral and bacterial concentrations

#Import libraries
#++++++++++++++++++++++++++++++++++++++++++++
import seaborn as sns
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
#++++++++++++++++++++++++++++++++++++++++++++


#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Lotka_Volterra(t,y):    
    return [r*y[0] - d*y[0]*y[1],
           c*d*y[0]*y[1] - m*y[1]]


def Lotka_Volterra_Logistic(t,y):    
    return [r*(1-(y[0]/K))*y[0] - d*y[0]*y[1],
           c*d*y[0]*y[1] - m*y[1]]

#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Min_Volume_B(t,y):

    return y[0] - 1
Min_Volume_B.terminal=True
Min_Volume_B.direction = 1
#:::::::::::::::::::::::::::::::::::::::::::::

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. MODEL AND SOLVER PARAMETERS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Parameters
#==============================
Dominant_Timescale='mr_Kick_Experiment_Logistic' #r/m/mr_kick_Experiment
Parameters,y0,time_0,time_f=Initial_Configuration(Dominant_Timescale)
#time_f_eq=0.1*time_f
time_f_eq=100
print(Parameters)
print(time_f)
print(y0)
r=Parameters['r'];d=Parameters['d'];c=Parameters['c']
m=Parameters['m'];K=Parameters['K']
y0=[m/(c*d), (r/d)*(1- m/(K*c*d))]
#==============================

#1.3. Solver Parameters
#==============================
Resolution=100;Step_Size=1/Resolution
time_eq=np.arange(time_0,time_f_eq,Step_Size)
step=1
#==============================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. MODEL SOLUTION
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.1. Initial conditions and Events
#==============================
Events_Model=[Min_Volume_B]
#==============================

#2.2. Solver Equilibrium
#====================================================================
sol_LV_EQ=solve_ivp(Lotka_Volterra,[time_eq[0],time_eq[-1]],y0,\
method='RK45',dense_output=True,events=Events_Model,max_step=step)

#2.3. Change time vector in case event triggered
#=====================================
try:
    final_t=sol_LV_EQ.t_events[0][0]
    time=time[:int(final_t)]
    print("Event found")
except IndexError:
    print('No event found')
#=====================================

#2.4. Calculate individual terms
#Discrete solution. Overlaps with time vector
z_eq = sol_LV_EQ.sol(time_eq)
Bacteria=z_eq[0][-1];Phage=z_eq[1][-1]
#====================================================================

#2.3. Solver Kick
#====================================================================

#Orders of magnitude
Orders_of_magnitude=np.logspace(-2,1,4,10)
Kick_list=[i*Phage for i in Orders_of_magnitude]

if Dominant_Timescale=='r_Kick_Experiment' or Dominant_Timescale=='r_Kick_Experiment_Logistic': 
    Kick={'Bacteria':[0.0]*len(Kick_list), 'Phage':Kick_list}
elif Dominant_Timescale=='m_Kick_Experiment' or Dominant_Timescale=='m_Kick_Experiment_Logistic':
#    Kick={'Bacteria':0.2, 'Phage':0.0}
    Kick={'Bacteria':Kick_list, 'Phage':[0.0]*len(Kick_list)}
elif Dominant_Timescale=='mr_Kick_Experiment' or Dominant_Timescale=='mr_Kick_Experiment_Logistic':
    Kick={'Bacteria':[0.0]*len(Kick_list), 'Phage':Kick_list}

time_K=np.arange(time_eq[-1],time_f,Step_Size)

Amplitudes_bacteria=[]
Amplitudes_phage=[]

for k in range(len(Kick_list)):
    print(Orders_of_magnitude[k])
    print(Kick['Bacteria'][k])
    print(Kick['Phage'][k])
    
    y_0_K=[Bacteria + Kick['Bacteria'][k]*Bacteria, Phage + Kick['Phage'][k]*Phage]
    print(y_0_K)

    sol_LV_K=solve_ivp(Lotka_Volterra,[time_K[0],time_K[-1]],y_0_K,\
method='RK45',dense_output=True,events=None,max_step=step)
    z_K = sol_LV_K.sol(time_K)
    z=np.concatenate((z_eq, z_K), axis=1)
    print(k)

    Min_bacteria=np.min(z[0])
    Max_bacteria=np.max(z[0])

    Min_phage=np.min(z[1])
    Max_phage=np.max(z[1])
    
    print(Min_bacteria, Max_bacteria)
    print(Min_phage, Max_phage)
    Amplitude_bacteria=Max_bacteria - Min_bacteria
    Amplitudes_bacteria.append(Amplitude_bacteria)
    
    Amplitude_phage=Max_phage - Min_phage
    Amplitudes_phage.append(Amplitude_phage)
    print('\n')
    
#y_0_K=[Bacteria + Kick['Bacteria']*Bacteria, Phage + Kick['Phage']*Phage]

plt.plot(Kick_list,Amplitudes_phage)
plt.show()

sol_LV_K=solve_ivp(Lotka_Volterra,[time_K[0],time_K[-1]],y_0_K,\
method='RK45',dense_output=True,events=None,max_step=step)
z_K = sol_LV_K.sol(time_K)

z=np.concatenate((z_eq, z_K), axis=1)
time=np.concatenate((time_eq,time_K))
#====================================================================

#2.4.3. Compute rates and extract individual terms
Bioprocesses=PerCapita_Rates_LV(z,time,Parameters,Dominant_Timescale)
Bioprocesses['Predation']=Bioprocesses.pop('Infection')
#====================================================================

#2.5. Convert Information to Dataframes
#===========================================
epsilon_1=0.1; epsilon_2=1
epsilon=epsilon_2

Bioprocesses_df=pd.DataFrame(data=Bioprocesses)
Total_Active_Bioprocesses,Initial_Dynamics=Get_Active_Bioprocesses(Bioprocesses_df,epsilon)

#print(Total_Active_Bioprocesses)
print(Initial_Dynamics)
#===========================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3.1. Configuration Figures
#=====================================================
#Save data
Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/Paper_Matt/'

Bioprocesses_df.to_csv=(Output_Path+'Data_colorbar.csv')
Bioprocesses_df_transp = Bioprocesses_df.transpose()

Max_val_0=Bioprocesses_df.max();Max_val=np.max(Bioprocesses_df.max())
Min_val_0=Bioprocesses_df.min();Min_val=np.min(Bioprocesses_df.min())

#Size of figures
cm = 1/2.54  # centimeters in inches
Width=8*cm;Height=5.0*cm #Width and height of plots

#Fontsizes
size_axis=7;size_ticks=6;size_title=5

#Gridspec parameters
Rows=1;Cols=2

#Extensions to save
Extensions=['.svg']

#Colors
cmap='RdBu';cmap_pieces=matplotlib.cm.get_cmap(cmap)
Phage_color=cmap_pieces(0.1);Bacteria_color=cmap_pieces(0.9)

#Linewidth
width_line=1
#=======================


#Plot Figures
#====================================================================
#Associate time vector with dataframe
Time_Bioprocesses_df=Bioprocesses_df.index.values.tolist()
Bioprocesses_df.set_index(time,inplace=True)

#Lotka volterra
#===============================================

#Figure
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
fig=figure(figsize=(Width, Height), dpi=300)

#Gridspec version
Gridspec_layout=False
if Gridspec_layout==True:
    gs=gridspec.GridSpec(Rows,Cols,width_ratios=[1,0.04])
    gs.update(left=0.1,right=0.7,bottom=0.08,top=0.97,wspace=0.1,hspace=0.1)
    ax_00=plt.subplot(gs[0,0])
    
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#2. Plot Dynamics
plt.plot(time, z[0].T,linewidth=width_line,color=Bacteria_color)
plt.plot(time, z[1].T,linewidth=width_line,color=Phage_color)

Min_bacteria=np.min(z[0])
Max_bacteria=np.max(z[0])

Min_phage=np.min(z[1])
Max_phage=np.max(z[1])

print(Min_bacteria, Max_bacteria)
print(Min_phage, Max_phage)

#Axes and legend
plt.xlim([time_0, time_f]);plt.ylim([1,5e8])
plt.xlim([1, time_f]);plt.ylim([1,5e8])
plt.ylabel('Concentration $\mathrm{(ml^{-1})}$',fontsize=size_axis)
plt.xlabel('Time (h)',fontsize=size_axis)
plt.legend(['Prey','Predator'],loc='best',fontsize=size_ticks,frameon=False)

#Ticks
step_ticks=int(time_f/10)
xticks_labels=[label for label in range(time_0,time_f+1,step_ticks)]
#xticks_labels=[label for label in np.arange(time_0,time_f+1,step_ticks)]

#plt.xticks(xticks_labels,fontsize=size_ticks)

plt.axvline(time_eq[-1], color='k', linewidth=0.75,linestyle='dotted')

plt.axhline(Max_bacteria, color='b', linewidth=0.5,linestyle='dotted')
plt.axhline(Min_bacteria, color='b', linewidth=0.5,linestyle='dotted')

plt.axhline(Max_phage, color='r', linewidth=0.5,linestyle='dotted')
plt.axhline(Min_phage, color='r', linewidth=0.5,linestyle='dotted')



    
#plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)
plt.xscale('log')
plt.yscale('log')

sns.despine(top=True, right=True, left=False, bottom=False) 

Name_Fig=str(Dominant_Timescale)
[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]

#plt.show()
#===============================================

'''
#Plot sum of mechanisms
#====================================================================
#Figure
fig=figure(figsize=(Width,Height),dpi=300)

#Optional gridspec version
Gridspec_layout=False
if Gridspec_layout==True:
    gs=gridspec.GridSpec(Rows,Cols,width_ratios=[1,0.04])
    gs.update(left=0.1,right=0.7,bottom=0.08,top=0.97,wspace=0.1,hspace=0.1)
    ax_00=plt.subplot(gs[0,0])

#Ticks
Scale='linear'
if Scale=='linear':
    step_ticks=int(time_f/10)
    xticks_labels=[label for label in range(time_0,time_f+1,step_ticks)] 
    plt.xticks(xticks_labels)
    plt.xticks(fontsize=size_ticks)
elif Scale=='log':
    plt.xscale('log')

plt.yticks(fontsize=size_ticks)
    
#Plot
plt.plot(time,Total_Active_Bioprocesses,color='k',linewidth=width_line)
plt.axvline(time_eq[-1], color='k', linewidth=1,linestyle='dotted') 

#Axes
plt.xlim([time_0,time_f]);plt.ylim([0,4])
plt.xlabel('Time (h)',fontsize=size_axis)
plt.ylabel('Number of active terms',fontsize=size_axis) 

sns.despine(top=True, right=True, left=False, bottom=False) 

Name_Fig='Sum_of_Mechanisms_'+str(Dominant_Timescale)
[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]

plt.show()


#Gridspec of heatmap and colbar
#===============================================
#1.Figure
fig=figure(figsize=(Width, Height),dpi=300)

#Gridspec version
Gridspec_layout=False
if Gridspec_layout==True:
    gs=gridspec.GridSpec(Rows,Cols,width_ratios=[1,0.04])
    gs.update(left=0.1,right=0.7,bottom=0.08,top=0.97,wspace=0.1,hspace=0.1)
    ax_00=plt.subplot(gs[0,0])
else:
    gs=gridspec.GridSpec(Rows,Cols,width_ratios=[1,0.00])
    gs.update(left=0.1,right=0.95,bottom=0.08,top=0.97,wspace=0.1,hspace=0.1)
    ax_00=plt.subplot(gs[0,0])

#3.4.1. Ticks, labels, and axes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Scale='linear'
if Scale=='linear':
    step_ticks=int(time_f/10)
    xticks_pos=[position for position in np.arange(Time_Bioprocesses_df[0],Time_Bioprocesses_df[-1]+1,step_ticks*Resolution)]
    xticks_labels=[label for label in np.arange(time_0,time_f+1,step_ticks)]

    #Reduce the number of columns in the dataframe and weight of svg
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Light_Heatmap_Vector=time[0::100]
    Light_Bioprocesses_df=Bioprocesses_df.loc[Light_Heatmap_Vector]
    Light_Bioprocesses_df_transp=Light_Bioprocesses_df.transpose()
    Light_time=Light_Bioprocesses_df.index.values.tolist()
    Light_time=np.asarray(Light_time)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    step_ticks=int(time_f/10)
    xticks_pos=[position for position in np.arange(Light_time[0],Light_time[-1]+1,step_ticks*Resolution/100)]
    print(xticks_pos)
    xticks_labels=[int(label) for label in np.arange(time_0,time_f+1,step_ticks)]
    print(xticks_labels)
    
    ax_heatmap=sns.heatmap(Light_Bioprocesses_df_transp,vmin=Min_val,vmax=Max_val,cmap="Reds",cbar=False)
    ax_00.set_xticks(xticks_pos)
    ax_00.set_xticklabels(xticks_labels,fontsize=size_ticks,rotation='horizontal')
    
elif Scale=='log':
    Logscale_Vector=Log_Sampling_Vector_fun(time,3,1)
    #Use log vector to logarithmically sample the dataset of the heatmap
    Log_Bioprocesses_df=Bioprocesses_df.loc[Logscale_Vector]
    Log_Bioprocesses_df_transp=Log_Bioprocesses_df.transpose()
    Log_Time_Bioprocesses_df=Log_Bioprocesses_df.index.values.tolist()
    
    Log_Time_Bioprocesses_df=np.asarray(Log_Time_Bioprocesses_df)
    xticks_pos=[Threshold_Finder(np.asarray(Log_Time_Bioprocesses_df), 0.1)[0],Threshold_Finder(Log_Time_Bioprocesses_df,1)[0], Threshold_Finder(Log_Time_Bioprocesses_df,10)[0], Threshold_Finder(Log_Time_Bioprocesses_df,100)[0]]

    xticks_labels=[0.1,1, 10,100]
    ax_heatmap=sns.heatmap(Log_Bioprocesses_df_transp,vmin=Min_val,vmax=Max_val,cmap="Reds",cbar=False)
    ax_00.set_xticks(xticks_pos)
    ax_00.set_xticklabels(xticks_labels,fontsize=size_ticks,rotation='horizontal')

    
#plt.xticks(rotation='horizontal',fontsize=size_ticks)
plt.yticks(fontsize=size_ticks)
ax_heatmap.tick_params(axis=u'y', which=u'both',length=0)
plt.xlabel('Time (h)',fontsize=size_axis)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
print(list(Kick.keys()))
Name_Fig='Heatmap_'+str(Dominant_Timescale)
[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]

plt.show()
#===============================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

'''
