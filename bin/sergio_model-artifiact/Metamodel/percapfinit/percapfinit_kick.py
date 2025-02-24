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


from scipy.signal import find_peaks
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

#Input
#=======================================
path='../../../data/input_output_data/'

file_name='Input_File_Kick.csv'
inputs=pd.read_csv(path+file_name)
inputs=inputs.set_index('description')
print(inputs)

scenario='decay_80_Beq_kick'

r=inputs.loc[scenario]['r']
d=inputs.loc[scenario]['a']
c=inputs.loc[scenario]['c']
m=inputs.loc[scenario]['m']
tau=inputs.loc[scenario]['tau']
time_f=inputs.loc[scenario]['t_f']

B_0_k=inputs.loc[scenario]['B0_kick']
P_0_k=inputs.loc[scenario]['P0_kick']
B_c=1/(c*d*tau);P_c=1/(d*tau);

time_f_eq=inputs.loc[scenario]['t_kick']
K=inputs.loc[scenario]['K']

if math.isnan(K)==True:
    print("hola")
    Beq=m/(c*d);Peq=(r/d)
else:
    Beq=m/(c*d);Peq=(r/d)*(1-m/(K*c*d))

y0=[Beq,Peq]
y_0_K=[B_0_k,P_0_k]
print("New_approach")
print(y_0_K)

B_c=1/(c*d*tau)
P_c=1/(d*tau) 
#==============================

#1.3. Solver Parameters
#==============================
time_0=0
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
if math.isnan(K)==True:
    sol_LV_EQ=solve_ivp(Lotka_Volterra,[time_eq[0],time_eq[-1]],y0,\
method='RK45',dense_output=True,events=Events_Model,max_step=step)
else:
    sol_LV_EQ=solve_ivp(Lotka_Volterra_Logistic,[time_eq[0],time_eq[-1]],y0,method='RK45',dense_output=True,events=Events_Model,max_step=step)
    
    
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
# print("Old_approach")
# print(y_0_K)

time_K=np.arange(time_eq[-1],time_f,Step_Size)

if math.isnan(K)==True:
    sol_LV_K=solve_ivp(Lotka_Volterra,[time_K[0],time_K[-1]],y_0_K,method='RK45',dense_output=True,events=None,max_step=step)
else:
    sol_LV_K=solve_ivp(Lotka_Volterra_Logistic,[time_K[0],time_K[-1]],y_0_K,method='RK45',dense_output=True,events=None,max_step=step)
    
z_K = sol_LV_K.sol(time_K)

z=np.concatenate((z_eq, z_K), axis=1)
time=np.concatenate((time_eq,time_K))

Solution_kick_dict={'time':time,'Bacteria':z[0],'Phage':z[1]}
Solution_kick_df=pd.DataFrame(data=Solution_kick_dict)

Save_csv=False
if Save_csv==True:
    Solution_kick_df.to_csv(path+'Solution_'+ scenario+'.csv')

#Solution_kick_df_slice=Solution_kick_df[(Solution_kick_df['time']>=time_K[0]) & (Solution_kick_df['time']<=time_K[0]+300)]

Solution_kick_df_slice=Solution_kick_df[(Solution_kick_df['time']>=time_K[0])]

Max_bacteria=Solution_kick_df_slice['Bacteria'].max()
Min_bacteria=Solution_kick_df_slice['Bacteria'].min()
Amplitude_bacteria=Max_bacteria - Min_bacteria

Minimum_bacteria_eq=Beq-Min_bacteria

Max_phage=Solution_kick_df_slice['Phage'].max()
Min_phage=Solution_kick_df_slice['Phage'].min()
Amplitude_phage=Max_phage - Min_phage

Minimum_phage_eq=Peq-Min_phage


#Find minima of bacteria
Negative_bacteria=-Solution_kick_df['Bacteria'] #find maxima of negative
peaks_bacteria, _ = find_peaks(Negative_bacteria, height=-1.01*Min_bacteria)
#Third minimum corresponds to double oscillation
first_min_B=peaks_bacteria[0]
time_first_min_B=Solution_kick_df_slice.loc[first_min_B]['time']
try:
    second_min_B=peaks_bacteria[1]
    time_second_min_B=Solution_kick_df_slice.loc[second_min_B]['time']
    period_B=time_second_min_B-time_first_min_B
    double_period_B=2*period_B
    relative_double_period_B=double_period_B/tau

    print("Double period bacteria")
    print(double_period_B)

    print("Relative period bacteria")
    print(relative_double_period_B)

except IndexError:
    print("no second minimum")



#Find minima of phage
Negative_phage=-Solution_kick_df['Phage'] #find maxima of negative
peaks_phage, _ = find_peaks(Negative_phage, height=-1.01*Min_phage)
#Third minimum corresponds to double oscillation
try:
    first_min_P=peaks_phage[0]
    time_first_min_P=Solution_kick_df_slice.loc[first_min_P]['time']
except IndexError:
    pass
try:
    second_min_P=peaks_phage[1]
    time_second_min_P=Solution_kick_df_slice.loc[second_min_P]['time']

    period_P=time_second_min_P-time_first_min_P
    double_period_P=2*period_P
    relative_double_period_P=double_period_P/tau

    print("Double period phage")
    print(double_period_P)

    print("Relative period phage")
    print(relative_double_period_P)
except IndexError:
    print("no second minimum phage")


print("\n")

print("Outputs of scenario " + str(scenario))
print("Final values")
print(Solution_kick_df.iloc[-1])
print("Amplitude bacteria")
print(Amplitude_bacteria)
print("Amplitude phage")
print(Amplitude_phage)

print("Beq-Bmin")
print(Minimum_bacteria_eq)
print("Peq-Pmin")
print(Minimum_phage_eq)

print("Maximum Bacteria")
print(Max_bacteria)

print("Minimum Bacteria")
print(Min_bacteria)

print("Maximum Phage")
print(Max_phage)
print("Minimum Phage")
print(Min_phage)

#====================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3.1. Configuration Figures
#=====================================================
#Save data
Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/Paper_Matt/'


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
width_line=2
#=======================

#Plot Figures
#====================================================================
#Associate time vector with dataframe
# Time_Bioprocesses_df=Bioprocesses_df.index.values.tolist()
# Bioprocesses_df.set_index(time,inplace=True)

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


#Axes and legend
plt.xlim([time_0, time_f]);plt.ylim([1,5e8])
plt.xlim([1, time_f]);
plt.ylabel('Concentration $\mathrm{(ml^{-1})}$',fontsize=size_axis)
plt.xlabel('Time (h)',fontsize=size_axis)
plt.legend(['Prey','Predator'],loc='best',fontsize=size_ticks,frameon=False)
plt.axvline(time_eq[-1], color='k', linewidth=0.75,linestyle='dotted')

plt.axhline(B_c, color=Bacteria_color, linewidth=0.75,linestyle='dotted')

plt.axhline(P_c, color=Phage_color, linewidth=0.75,linestyle='dotted')

#Ticks
step_ticks=int(time_f/10)
#xticks_labels=[label for label in range(time_0,time_f+1,step_ticks)]
xticks_labels=[label for label in np.arange(time_0,time_f+1,step_ticks)]

#plt.xticks(xticks_labels,fontsize=size_ticks)

#plt.xticks(fontsize=size_ticks)
if time_f>1e4:
    plt.xscale('log')
    plt.xlim([1, time_f]);
    plt.xticks(fontsize=size_ticks)
else:
    plt.xticks(xticks_labels,fontsize=size_ticks)
    
plt.yscale('log')    
plt.yticks(fontsize=size_ticks)

    


sns.despine(top=True, right=True, left=False, bottom=False) 

Name_Fig=str(scenario)
[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]


plt.show()
#===============================================

'''
=======
plt.show()
#===============================================

>>>>>>> 9944fffb8a9fa6512692733ae4d9ae95bfad9942
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
