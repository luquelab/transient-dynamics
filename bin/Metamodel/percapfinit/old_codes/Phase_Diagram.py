#7/3/2023. This code plot the phase diagram of the results of the metamodel analysis.

#Import libraries
#++++++++++++++++++++++++++++++++++++++++++++
import seaborn as sns
import numpy as np
from scipy.integrate import odeint
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.pyplot as plt
from decimal import Decimal
import math
from scipy.integrate import solve_ivp
import pandas as pd
from matplotlib.pyplot import figure
import copy
from Metamodel_Functions import *
import sys
#++++++++++++++++++++++++++++++++++++++++++++

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#0. Lotka-Volterra Equations
#=============================================================
def Lotka_Volterra(t,y):    

    return [r*y[0] - d*y[0]*y[1],
           c*d*y[0]*y[1] - m*y[1]]

#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Min_Volume_B(t,y):

    return y[0] - 1
Min_Volume_B.terminal=True
Min_Volume_B.direction = 1
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================

#m dominant scale
#====================================================================
#1. Burst and decay active
#=============================================================
def Lotka_Volterra_BD(t,y,Thresholds):    
    return [0,
           c*d*y[0]*y[1] - m*y[1]]

#Activate infection
def Switch_Predation_On(t,y,Thresholds):

    return y[1] - Thresholds
Switch_Predation_On.terminal=True
Switch_Predation_On.direction = 0

def Switch_Burst_On(t,y,Thresholds):
    return y[0] - Thresholds

Switch_Burst_On.terminal=True
Switch_Burst_On.direction = 0
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================

#2. Predation, burst and decay active
#=============================================================
def Lotka_Volterra_PBD(t,y,Thresholds):    
    return [-d*y[0]*y[1],
           c*d*y[0]*y[1] - m*y[1]]

#Inactivation of burst
def Switch_Burst_Off(t,y,Thresholds):
    return y[0] - Thresholds
Switch_Burst_Off.terminal=True
Switch_Burst_Off.direction = 0
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================

#m and r dominant scales
#====================================================================

#1. Growth, decay, and burst active
#=============================================================
def Lotka_Volterra_GBD(t,y,Thresholds):    
    return [r*y[0] ,
           c*d*y[0]*y[1] - m*y[1]]
#=============================================================

#2. Growth, decay, burst, and predation active
#=============================================================
def Lotka_Volterra_GPBD(t,y,Thresholds):
    return [r*y[0] - d*y[0]*y[1],
           c*d*y[0]*y[1] - m*y[1]]
#=============================================================

#3. Growth, decay, and predation active
#=============================================================
def Lotka_Volterra_GPD(t,y,Thresholds):    
    return [r*y[0] - d*y[0]*y[1],
            - m*y[1]]
#=============================================================

#4. Growth and decay
#=============================================================
def Lotka_Volterra_GD(t,y,Thresholds):    
    return [r*y[0],-m*y[1]]
#=============================================================

#5. Growth
#=============================================================
def Lotka_Volterra_G(t,y,Thresholds):    
    return [r*y[0],0]
#=============================================================

#6. Growth and Burst
#=============================================================
def Lotka_Volterra_GB(t,y,Thresholds):    
    return [r*y[0],c*d*y[0]*y[1]]
#=============================================================

#7. Growth,Predation and Burst
#=============================================================
def Lotka_Volterra_GPB(t,y,Thresholds):    
    return [r*y[0]-d*y[0]*y[1],c*d*y[0]*y[1]]
#=============================================================

#8. Growth, and Predation
#=============================================================
def Lotka_Volterra_GP(t,y,Thresholds):    
    return [r*y[0]-d*y[0]*y[1],0]
#=============================================================

#9. Predation, Burst, and Decay
#=============================================================
def Lotka_Volterra_PBD(t,y,Thresholds):    
    return [-d*y[0]*y[1],
           c*d*y[0]*y[1] - m*y[1]]
#=============================================================

#10. Predation and Decay
#=============================================================
def Lotka_Volterra_PD(t,y,Thresholds):    
    return [-d*y[0]*y[1],
            - m*y[1]]
#=============================================================

#11. Decay
#=============================================================
def Lotka_Volterra_D(t,y,Thresholds):    
    return [0,- m*y[1]]
#=============================================================

#Events
#=======================================
#Switch Burst
def Switch_Burst(t,y,Thresholds):
    return y[0] - Thresholds
Switch_Burst.terminal=True
Switch_Burst.direction = 0

#Switch Predation
def Switch_Predation(t,y,Thresholds):
    return y[1] - Thresholds
Switch_Predation.terminal=True
Switch_Predation.direction = 0
#=======================================

#5. Time vector for simplified dynamics
#=============================================================
def Build_Time_Vector(t_0, t_f):

    Resolution=100;Steps=time_f*Resolution + 1
    Step_Size=1/Resolution
    Time_Vector=np.arange(t_0,t_f,Step_Size)

    return Time_Vector
#=============================================================

#6. Solution for simplified dynamics
#=============================================================
def Solve_Simplified_Model(y0_fun0,t0_fun,tf_fun,Event_thresholds,Model_fun,Event_functions):

    #Build time vector and set initial conditions
    step=0.01    
    t_in=Build_Time_Vector(t0_fun,tf_fun)
    y0_fun=[ y0_fun0[0], y0_fun0[1] ]

    #Solve simplified model
    sol_Model_Raw=\
    solve_ivp(lambda t,x:Model_fun(t,x,Event_thresholds),\
    [t_in[0],t_in[-1]],\
    y0_fun,method='RK45',\
    dense_output=True,\
    events=[lambda t,x:Event_functions[0](t,x,Event_thresholds[0]),\
            lambda t,x:Event_functions[1](t,x,Event_thresholds[1])],
              max_step=step)

    #Manage events
    #------------------------------------------------------
    Active_Event_Types=[0]*len(Event_functions) 

    #Remove event if it is equal to starting time
    #-------------------------------------------------------
    for index in range(len(Event_functions)):
        if t_in[0] in sol_Model_Raw.t_events[index]:
            sol_Model_Raw.t_events[index]=\
            np.delete(sol_Model_Raw.t_events[index],0)
            sol_Model_Raw.y_events[index]=\
            np.delete(sol_Model_Raw.y_events[index],0,0)
    #-------------------------------------------------------

    #Count number of active events
    #-------------------------------------------------------
    for item in range(len(sol_Model_Raw.t_events)):
        if len(sol_Model_Raw.t_events[item])>0:
            Active_Event_Types[item]=1

    Types_of_Events=sum(Active_Event_Types)
    #-------------------------------------------------------

#Move on if no events were found. Find earliest event otherwise
    if Types_of_Events==0:
        print("No event found")
        t_out=t_in
        sol_Model=sol_Model_Raw.sol(t_out)
        y_f=[sol_Model[0][-1], sol_Model[1][-1]]
        ind_first_event=None
            
    else:
        #Take the earliest times of events
        #-----------------------------------------------
        Earliest_events={}
        for index in range(len(sol_Model_Raw.t_events)):
#            print(index)
            try:
                Earliest_events[index]=sol_Model_Raw.t_events[index][0]
            except IndexError:
                continue

        print(Earliest_events)
        ind_first_event=min(Earliest_events, key=Earliest_events.get)
        #-----------------------------------------------

        #Declare final time, final conditions and new time vector
        #------------------------------------------------
        t_f=sol_Model_Raw.t_events[ind_first_event][0]
        y_f=sol_Model_Raw.y_events[ind_first_event][0]
        t_out=Build_Time_Vector(t_in[0],t_f)
        #------------------------------------------------
        
        #Solution of simplified dynamics
        #----------------------------------
        sol_Model=sol_Model_Raw.sol(t_out)
        #----------------------------------
        
        print('Event found')
        print(t_f)
        print(y_f)
                        
    return sol_Model, t_out, y_f, ind_first_event
#=============================================================

#7. Monitor active/inactive dynamics
#=============================================================
def Monitor_Dynamics(Dynamics_Dict,Index_var):

    if Index_var==0:
        Dynamics_Dict['Predation']=Dynamics_Dict['Predation']-1
        Dynamics_Dict['Predation']=abs(Dynamics_Dict['Predation'])

    elif Index_var==1:
        Dynamics_Dict['Burst']=Dynamics_Dict['Burst']-1
        Dynamics_Dict['Burst']=abs(Dynamics_Dict['Burst'])  
    
    else:        
        Dynamics_Dict=Dynamics_Dict
    
    Hash_Dynamics_0=[]    
    for term in Terms_Ordered:                                    
        Hash_Dynamics_0.append(Dynamics_Dict[term])
        Hash_Dynamics=''.join(str(value) for value in Hash_Dynamics_0)

    
    return Hash_Dynamics, Dynamics_Dict
#=============================================================

#Write output to file
#=============================================================
def write_to_file(File, Name):
    
    with open(Name,'w', encoding="utf-8") as file_to_save:
        for item in zip(File[0],File[1]):
            file_to_save.write( str(item) + "\n")
        
    return None
#=============================================================

#=============================================================
def Get_hash_of_initial_dynamics(Dynamics,Terms_Ordered_fun):

    Active_Terms_List=[]

    for term in Terms_Ordered_fun:
        Active_Terms_List.append(Dynamics[term])

    Initial_Hash=''.join(str(value) for value in Active_Terms_List)

    return Initial_Hash
#=============================================================

#Concatenated simplified dynamics
#=============================================================
def Concatenated_simplified_dynamics(t0_fun, tf_fun, y0_fun, step_fun, Initial_dynamics_fun, tipping_points ,Dictionary_dynamics_fun ,hash_dynamics_fun,list_events):

    Simplified_dynamics_list=[];Simplified_times_list=[]
    counter=0

    while t0_fun<tf_fun-step_fun:
        
        #Solve individual dynamic
        #---------------------------------------------------------
        z_fun,time_fun,y_f_fun,type_event=Solve_Simplified_Model(y0_fun,t0_fun,tf_fun,tipping_points,Dictionary_dynamics_fun[hash_dynamics_fun],[list_events[0], list_events[1]])
        #-----------------------------------------------------------
        #Update information for next iterations
        #-----------------------------------------------------------
        t0_fun=time_fun[-1];y0_fun=y_f_fun;counter+=1
        hash_dynamics_fun,Initial_dynamics_fun=Monitor_Dynamics(Initial_dynamics_fun,type_event)
        #-----------------------------------------------------------

        #Store information of this iteration
        #--------------------------------------
        Simplified_dynamics_list.append(z_fun)
        Simplified_times_list.append(time_fun)
        #--------------------------------------

        #------------------------------------------------------------
    Simplified_dynamics_fun=np.concatenate(Simplified_dynamics_list,axis=1)
    Simplified_times_fun=np.concatenate(Simplified_times_list)
    Critical_Times=[time[-1] for time in Simplified_times_list]
    #.............................................................

    return Simplified_dynamics_fun, Simplified_times_fun, Critical_Times
#=============================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. MODEL AND SOLVER PARAMETERS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Choose dominant timescale with paramater values associated
#====================================================================
Dominant_Timescale='m' #r, m, mr_Equ, mr_Disequ
Parameters,y0,time_0,time_f=Initial_Configuration(Dominant_Timescale)
print(Parameters)
print(time_f)
print(y0)
r=Parameters['r'];d=Parameters['d'];c=Parameters['c']
m=Parameters['m'];K=Parameters['K']
#====================================================================

#1.3. Build time vector
#==============================
Resolution=100;Steps=time_f*Resolution + 1
time=np.linspace(time_0,time_f,Steps)
step=0.01
#==============================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. MODEL SOLUTION
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Solve full Model with event
#===================================================================
sol_LV=solve_ivp(Lotka_Volterra,[time[0],time[-1]],y0,\
method='RK45',dense_output=True,events=Min_Volume_B,max_step=step)

try:
    final_t=sol_LV.t_events[0][0]
    time=time[:int(final_t)]
except IndexError:
    print('No event found')
#===================================================================

#2.4. Calculate individual terms
#====================================================================
#2.4.1. Discrete solution. Overlaps with time vector
z = sol_LV.sol(time)

#Find total active processes, critic concentrations, and times
Bioprocesses=PerCapita_Rates_LV(z,time,Parameters,Dominant_Timescale)
Bioprocesses['Predation']=Bioprocesses.pop('Infection')
#Convert to dataframe
Bioprocesses_df=pd.DataFrame(data=Bioprocesses)
#Associate time vector with dataframe
Time_Bioprocesses_df=Bioprocesses_df.index.values.tolist()
#....................................................................

#Universal order of terms
Terms_Ordered=['Growth', 'Predation', 'Burst', 'Decay']

Dictionary_Dynamics={'1111':Lotka_Volterra_GPBD, '1011':Lotka_Volterra_GBD, '1101': Lotka_Volterra_GPD, '1001':Lotka_Volterra_GD, '1000':Lotka_Volterra_G, '1010':Lotka_Volterra_GB, '1100':Lotka_Volterra_GP, '1110':Lotka_Volterra_GPB, '0011':Lotka_Volterra_BD, '0111':Lotka_Volterra_PBD, '0101':Lotka_Volterra_PD, '0001': Lotka_Volterra_D}

#Define epsilon
epsilon_1=0.1;epsilon_2=1
epsilon=epsilon_2

epsilon_list=[0.1,0.5,1]

Simplified_dynamics_dict={}
Simplified_times_dict={}
Critical_times_dict={}
Critical_Concentrations={}
Active_Bioprocesses_dict={}

for epsilon in epsilon_list:
#Find active bioprocesses and initial dynamics
    Total_Active_Bioprocesses,Initial_Dynamics=Get_Active_Bioprocesses(Bioprocesses_df,epsilon)


    Active_Bioprocesses_dict[epsilon]=Total_Active_Bioprocesses
    Hash_Initial_Dynamics=Get_hash_of_initial_dynamics(Initial_Dynamics,Terms_Ordered)


    if Dominant_Timescale=='r':
        Critical_Phage=epsilon*r/d;Critical_Bact=epsilon*r/(c*d)

    elif Dominant_Timescale=='m':
        Critical_Phage=epsilon*m/d;Critical_Bact=epsilon*m/(c*d)

    elif Dominant_Timescale=='mr_Equ' or Dominant_Timescale=='mr_Disequ':
        Critical_Phage=epsilon*r/d;Critical_Bact=epsilon*r/(c*d)
    

    Initial_Dynamics_copy=Initial_Dynamics.copy()
    t0=time_0;y0_new=y0
    
    Crit_Concs=[Critical_Phage,Critical_Bact]
    Critical_Concentrations[epsilon]=Crit_Concs
        
    Simplified_dynamics, Simplified_times, Critical_times=Concatenated_simplified_dynamics(t0,time_f,y0_new,step,Initial_Dynamics_copy,Crit_Concs,Dictionary_Dynamics,Hash_Initial_Dynamics,[Switch_Predation,Switch_Burst])
        
    Simplified_dynamics_dict[epsilon]=Simplified_dynamics
    Simplified_times_dict[epsilon]=Simplified_times
    Critical_times_dict[epsilon]=Critical_times
    
    Rel_Error,Rel_Err_Bacteria,Rel_Err_Phage=Get_Relative_Error(z,Simplified_dynamics)
    
    print("Error")
    print(np.mean(Rel_Error), np.mean(Rel_Err_Bacteria), np.mean(Rel_Err_Phage))
#--------------------------------------------------------------------


#Plots
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#3.1. Configuration of plots
#====================================================================
#Path to save figure
Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/Paper_Matt/'

#Fontsizes
size_axis=7;size_ticks=6;size_title=5

#Figure Size
cm = 1/2.54  # centimeters in inches
Width=8*cm;Height=4*cm #Width and height of plots

#Gridspec parameters
Rows=1;Cols=2

#Colors
cmap='RdBu';cmap_pieces= matplotlib.cm.get_cmap(cmap)
Phage_color=cmap_pieces(0.1);Bacteria_color=cmap_pieces(0.9)

#Extensions to save
Extensions=['.png','.svg']

#Linewidth
width_line=1

#Line dashes
dash_lengths=[1,3,6]
dash_spaces=[1,2,4]

Line_type={}
counter=0
for eps in epsilon_list:
    Line_type[eps]=(dash_lengths[counter],dash_spaces[counter])
    counter+=1
#====================================================================


#Plot Full vs Simplified dynamics
#====================================================================
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


#Plot phase diagram
#====================================================================
plt.plot(z[1].T,z[0].T,color='black',linewidth=1,label='Full')

#count_dashes=0
for epsilon in epsilon_list:

    
    # Bacterial_oscillation_cycles=Threshold_Finder(Simplified_dynamics_dict[epsilon][0],y0[0]-1e-2)
    # Phage_oscillation_cycles=Threshold_Finder(Simplified_dynamics_dict[epsilon][1],y0[1])
    # First_cycle=[Bacterial_oscillation_cycles[1],Phage_oscillation_cycles[1]]
    # Longest_cycle=np.max(First_cycle)
    
    # Simplified_Bacteria_Cycle=Simplified_dynamics_dict[epsilon][0].T[:Longest_cycle]
    # Simplified_Phage_Cycle=Simplified_dynamics_dict[epsilon][1].T[:Longest_cycle]

    Max_value_phage=np.max(Simplified_dynamics_dict[epsilon][1].T)
    Min_value_phage=np.min(Simplified_dynamics_dict[epsilon][1].T)
    
    Max_value_bact=np.max(Simplified_dynamics_dict[epsilon][0].T)
    Min_value_bact=np.min(Simplified_dynamics_dict[epsilon][0].T)

    plt.hlines(Critical_Concentrations[epsilon][1],xmin=0.9*Min_value_phage,xmax=1.1*Max_value_phage,color=Bacteria_color,linewidth=1,linestyle='-')
    plt.vlines(Critical_Concentrations[epsilon][0],ymin=0.9*Min_value_bact,ymax=0.9*Max_value_bact,color=Phage_color,linewidth=1,linestyle='-')
    

#    plt.plot(Simplified_Phage_Cycle,Simplified_Bacteria_Cycle,color='gray',linestyle='dashed',dashes=Line_type[epsilon],linewidth=1,label='$\epsilon=$'+str(epsilon))

    plt.plot(Simplified_dynamics_dict[epsilon][1], Simplified_dynamics_dict[epsilon][0],color='gray',linestyle='dashed',dashes=Line_type[epsilon],linewidth=1,label='$\epsilon=$'+str(epsilon))


plt.xlabel('Predator concentration $\mathrm{(ml^{-1})}$',fontsize=size_axis)
plt.ylabel('Prey concentration $\mathrm{(ml^{-1})}$',fontsize=size_axis)


#Ticks
#step_ticks=10
#xticks_labels=[label for label in range(time_0,time_f+1,step_ticks)]
#plt.xticks(xticks_labels)

plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)
plt.ticklabel_format(axis='both', style='sci', scilimits=(4,6))

plt.legend(loc='upper right',fontsize=size_ticks,frameon=False)

sns.despine(top=True, right=True, left=False, bottom=False)

Name_Figure='Phase_Diagrams_'+str(Dominant_Timescale)

for ext in Extensions:
    plt.savefig(Output_Path+Name_Figure+ext,dpi=300)

plt.show()
#====================================================================

#Plot Full Dynamics
#====================================================================

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

plt.plot(time,z[0].T,color=Bacteria_color,linewidth=1,label='Prey')
plt.plot(time,z[1].T,color=Phage_color,linewidth=1,label='Predator') 


for epsilon in epsilon_list:
    plt.plot(Simplified_times_dict[epsilon],Simplified_dynamics_dict[epsilon][0].T,color='gray',linestyle='dashed',dashes=Line_type[epsilon],linewidth=1,label='$\epsilon=$'+ str(epsilon))

    plt.plot(Simplified_times_dict[epsilon],Simplified_dynamics_dict[epsilon][1].T,color='gray',linestyle='dashed',dashes=Line_type[epsilon],linewidth=1)

     
#Axes, title and ticks
plt.ylabel('Concentration $\mathrm{(ml^{-1})}$',fontsize=size_axis)
plt.xlabel('Time (h)',fontsize=size_axis)
plt.legend(loc='best',fontsize=size_ticks,frameon=False)

#Ticks
step_ticks=2
xticks_labels=[label for label in range(time_0,time_f+1,step_ticks)]
plt.xticks(xticks_labels)
plt.xticks(fontsize=size_ticks)
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks) 

sns.despine(top=True, right=True, left=False, bottom=False)

#x,y limits
plt.xlim([time[0],time[-1]])
plt.yscale('log')
plt.ylim([1, 10*np.max(z)])

if Dominant_Timescale=='m':
    plt.xscale('linear')

#Save and show plots
Name_Figure='Dynamics_eps_'+str(epsilon)+'_'+str(Dominant_Timescale)

for ext in Extensions:
    plt.savefig(Output_Path+Name_Figure+ext,dpi=300)
plt.show()
#====================================================================


#Plot sum of mechanisms
#====================================================================
fig=figure(figsize=(Width,Height),dpi=300)

#Gridspec version of figure
Gridspec_layout=False
if Gridspec_layout==True:
    gs=gridspec.GridSpec(Rows,Cols,width_ratios=[1,0.04])
    gs.update(left=0.1,right=0.7,bottom=0.08,top=0.97,wspace=0.1,hspace=0.1)
    ax_00=plt.subplot(gs[0,0])

#Plot
for epsilon in Active_Bioprocesses_dict:
    plt.plot(time,Active_Bioprocesses_dict[epsilon],color='k',linewidth=width_line,linestyle='--',dashes=Line_type[epsilon],label='$\epsilon=$'+ str(epsilon))
    
#Axes
plt.xlim([time_0,time_f]);plt.ylim([0,4.5])
plt.ylabel('Number of active terms',fontsize=size_axis)
plt.xlabel('Time (h)',fontsize=size_axis)
plt.legend(loc='best',fontsize=size_ticks,frameon=False)

#Ticks 
if Dominant_Timescale=='r':
    step_ticks=2
    xticks_labels=[label for label in range(time_0,time_f+1,step_ticks)]
    plt.xticks(xticks_labels)
    plt.xticks(fontsize=size_ticks)
elif Dominant_Timescale=='r':
    step_ticks=20
    xticks_labels=[label for label in range(time_0,time_f+1,step_ticks)]
    plt.xticks(xticks_labels)
    plt.xticks(fontsize=size_ticks) 
elif Dominant_Timescale=='m':
    plt.xscale('linear')
#    plt.xscale('log')

plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)
sns.despine(top=True, right=True, left=False, bottom=False)
Name_Fig='Sum_of_Mechs'+str(Dominant_Timescale)

[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]
plt.show()
#====================================================================

#Heatmap
#====================================================================
fig=figure(figsize=(Width, Height),dpi=300)

Bioprocesses_df_transp = Bioprocesses_df.transpose()
Max_val_0=Bioprocesses_df.max();Max_val=np.max(Bioprocesses_df.max())
Min_val_0=Bioprocesses_df.min();Min_val=np.min(Bioprocesses_df.min())

#Gridspec version of figure
Gridspec_layout=True
if Gridspec_layout==True:
    gs=gridspec.GridSpec(Rows,Cols,width_ratios=[1,0.04])
    gs.update(left=0.1,right=0.7,bottom=0.08,top=0.97,wspace=0.1,hspace=0.1)
    ax_00=plt.subplot(gs[0,0])

#3.4.1. Ticks, labels, and axes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
if Dominant_Timescale=='r' or Dominant_Timescale=='mr_Disequ':
    step_ticks=2 
    xticks_pos=[position for position in np.arange(Time_Bioprocesses_df[0],Time_Bioprocesses_df[-1]+1,step_ticks*Resolution)]
    xticks_labels=[label for label in np.arange(time_0,time_f+1,step_ticks)] 
    ax_heatmap=sns.heatmap(Bioprocesses_df_transp,vmin=Min_val,vmax=Max_val,cmap="Reds",cbar=False)
    ax_00.set_xticks(xticks_pos)
    ax_00.set_xticklabels(xticks_labels,fontsize=size_ticks,rotation='horizontal')

elif Dominant_Timescale=='m':
    #Linear sampling of the dataframe
    #--------------------------------------------------------------
    Light_Heatmap_Vector=time[0::100]
    Light_Bioprocesses_df=Bioprocesses_df.loc[Light_Heatmap_Vector]
    Light_Bioprocesses_df_transp=Light_Bioprocesses_df.transpose()
    Light_time=Light_Bioprocesses_df.index.values.tolist()
    Light_time=np.asarray(Light_time)

    step_ticks=int(time_f/10) 
    xticks_pos=[position for position in np.arange(Light_time[0],Light_time[-1]+1,step_ticks*Resolution/100)]
    print(xticks_pos) 
    xticks_labels=[int(label) for label in np.arange(time_0,time_f+1,step_ticks)]
    print(xticks_labels)

    ax_heatmap=sns.heatmap(Light_Bioprocesses_df_transp,vmin=Min_val,vmax=Max_val,cmap="Reds",cbar=False)
    #--------------------------------------------------------------
    ax_00.set_xticks(xticks_pos)
    ax_00.set_xticklabels(xticks_labels,fontsize=size_ticks,rotation='horizontal')
    
elif Dominant_Timescale=='mr_Equ' or Dominant_Timescale=='mr_Disequ': 
    ax_heatmap=sns.heatmap(Bioprocesses_df_transp,vmin=Min_val,vmax=Max_val,cmap="Reds",cbar=False)


#plt.xticks(rotation='horizontal',fontsize=size_ticks)
plt.yticks(fontsize=size_ticks)
ax_heatmap.tick_params(axis=u'y', which=u'both',length=0)
plt.xlabel('Time (h)',fontsize=size_axis)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Name_Fig='Terms_'+str(Dominant_Timescale)
[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]
plt.show()
#====================================================================

