#7/3/2023. This code plots the phase diagram of the results of the metamodel analysis.

#Import libraries
#++++++++++++++++++++++++++++++++++++++++++++
import seaborn as sns
import numpy as np
import scipy
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
from matplotlib import colormaps
#++++++++++++++++++++++++++++++++++++++++++++


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


#r on, m off
#====================================================================

#1. Growth
#=============================================================
def Lotka_Volterra_G(t,y,Thresholds):    
    return [r*y[0],0]
#=============================================================

#2. Growth and Burst
#=============================================================
def Lotka_Volterra_GB(t,y,Thresholds):    
    return [r*y[0],c*d*y[0]*y[1]]
#=============================================================

#3. Growth,Predation and Burst
#=============================================================
def Lotka_Volterra_GPB(t,y,Thresholds):    
    return [r*y[0]-d*y[0]*y[1],c*d*y[0]*y[1]]
#=============================================================

#4. Growth, and Predation
#=============================================================
def Lotka_Volterra_GP(t,y,Thresholds):    
    return [r*y[0]-d*y[0]*y[1],0]
#=============================================================

#====================================================================

#m on, r off
#====================================================================
#5. Burst and decay active
#=============================================================
def Lotka_Volterra_BD(t,y,Thresholds):    
    return [0,
           c*d*y[0]*y[1] - m*y[1]]

#Activate infection
def Switch_Predation_On(t,y,Thresholds):

    return y[1] - Thresholds

def Switch_Burst_On(t,y,Thresholds):
    return y[0] - Thresholds
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================

#6. Predation, Burst, and Decay
#=============================================================
def Lotka_Volterra_PBD(t,y,Thresholds):    
    return [-d*y[0]*y[1],
           c*d*y[0]*y[1] - m*y[1]]
#=============================================================

#7. Predation and Decay
#=============================================================
def Lotka_Volterra_PD(t,y,Thresholds):    
    return [-d*y[0]*y[1],
            - m*y[1]]
#=============================================================

#8. Decay
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
#====================================================================

#m on, and r on
#====================================================================

#9. Growth, decay, and burst active
#=============================================================
def Lotka_Volterra_GBD(t,y,Thresholds):    
    return [r*y[0] ,
           c*d*y[0]*y[1] - m*y[1]]
#=============================================================

#10. Growth, decay, burst, and predation active
#=============================================================
def Lotka_Volterra_GPBD(t,y,Thresholds):
    return [r*y[0] - d*y[0]*y[1],
           c*d*y[0]*y[1] - m*y[1]]
#=============================================================

#11. Growth, decay, and predation active
#=============================================================
def Lotka_Volterra_GPD(t,y,Thresholds):    
    return [r*y[0] - d*y[0]*y[1],
            - m*y[1]]
#=============================================================

#12. Growth and decay
#=============================================================
def Lotka_Volterra_GD(t,y,Thresholds):    
    return [r*y[0],-m*y[1]]
#=============================================================
#====================================================================

#m off, and r off
#====================================================================
#13. Burst and predation active
#=============================================================
def Lotka_Volterra_PB(t,y,Thresholds):
    return [- d*y[0]*y[1],
           c*d*y[0]*y[1]]
#=============================================================

#14. Predation active
#=============================================================
def Lotka_Volterra_P(t,y,Thresholds):
    return [- d*y[0]*y[1],0]
#=============================================================

#15. Burst active
#=============================================================
def Lotka_Volterra_B(t,y,Thresholds):
    return [0,c*d*y[0]*y[1]]
#=============================================================

#16. Nothing active
#=============================================================
def Lotka_Volterra_Not(t,y,Thresholds):
    return [0,0]
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

#====================================================================

#====================================================================

#Write output to file
#=============================================================
def write_to_file(File, Name):
    
    with open(Name,'w', encoding="utf-8") as file_to_save:
        for item in zip(File[0],File[1]):
            file_to_save.write( str(item) + "\n")
        
    return None
#=============================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. MODEL AND SOLVER PARAMETERS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Choose dominant timescale with paramater values associated
#====================================================================
path='../../../data/input_output_data/'
file_name='results_model_full.csv'
inputs=pd.read_csv(path+file_name)                              
inputs=inputs.set_index('description')                               
print(inputs)                                    

#scenarios: growth, decay, growth-decay,growth-decay_regimes, no-growth-no-decay
scenario='growth'

r=inputs.loc[scenario]['r']
d=inputs.loc[scenario]['a']                                          
c=inputs.loc[scenario]['c']                                          
m=inputs.loc[scenario]['m']
K=inputs.loc[scenario]['K']

Parameters={'r':r,'d':d,'c':c,'m':m,'K':K}

tau=inputs.loc[scenario]['tau']                                      
time_f=inputs.loc[scenario]['t_f']
time_0=inputs.loc[scenario]['t_0']
Beq=m/(c*d);Peq=r/d

B0=inputs.loc[scenario]['B0']
P0=inputs.loc[scenario]['P0']

y0=[B0,P0]

Active_order=1
#====================================================================

#1.3. Build time vector
#==============================
Resolution=100;Steps=time_f*Resolution + 1
time=np.linspace(time_0,time_f,Steps)
step=0.01
#==============================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.FULL MODEL SOLUTION
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Solve full Model with event
#===================================================================
sol_LV=solve_ivp(Lotka_Volterra,[time[0],time[-1]],y0,method='RK45',dense_output=True,events=Min_Volume_B,max_step=step)

try:
    final_t=sol_LV.t_events[0][0]
    time=time[:int(final_t)]
except IndexError:
    print('No event found')
#===================================================================

#Solution of the full model and per capita terms
#====================================================================
z = sol_LV.sol(time)

Solution_dict={'time': time, 'B_full':z[0].T , 'P_full': z[1].T }
Solution_df=pd.DataFrame(data=Solution_dict)

#Find where agents have a concentration of one per ml and save them
#..............................................................................................
One_per_ml_Indices,One_per_ml_concentrations=Find_one_per_ml(Solution_df, 'B_full', 'P_full',1)

One_per_ml_dict={}
for agent in One_per_ml_Indices:
    print(agent)
    for index in One_per_ml_Indices[agent]:
        print(Solution_df.iloc[index])
        print(Solution_df.iloc[index]['time'])
        try:
            One_per_ml_dict[agent].append(Solution_df.iloc[index]['time'])
        except KeyError:
            One_per_ml_dict[agent]=[Solution_df.iloc[index]['time']]  
#....................................................................


#Find critical concentrations, total active processes, and times
#................................................................
Critical_Phage=Active_order/(tau*d);Critical_Bact=Active_order/(c*d*tau)
Crit_Concs=[Critical_Phage,Critical_Bact]

Bioprocesses=PerCapita_Rates_Timescale(z,time,Parameters,time_f)
Bioprocesses['Predation']=Bioprocesses.pop('Infection')
Bioprocesses_df=pd.DataFrame(data=Bioprocesses)
#................................................................

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3. SIMPLIFIED MODEL SOLUTION
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Universal order of terms and dictionary of dynamic regimes
#................................................................
Terms_Ordered=['Growth', 'Predation', 'Burst', 'Decay']

Dictionary_Dynamics={'1111':Lotka_Volterra_GPBD, '1011':Lotka_Volterra_GBD, '1101': Lotka_Volterra_GPD, '1001':Lotka_Volterra_GD, '1000':Lotka_Volterra_G, '1010':Lotka_Volterra_GB, '1100':Lotka_Volterra_GP, '1110':Lotka_Volterra_GPB, '0011':Lotka_Volterra_BD, '0111':Lotka_Volterra_PBD, '0101':Lotka_Volterra_PD, '0001': Lotka_Volterra_D, '0010':Lotka_Volterra_B, '0110':Lotka_Volterra_PB, '0100':Lotka_Volterra_P, '0000':Lotka_Volterra_Not}
#................................................................

#Find number of active bioprocesses and initial dynamics
#................................................................
Total_Active_Bioprocesses,Initial_Dynamics=Get_Active_Bioprocesses(Bioprocesses_df,Active_order)
Hash_Initial_Dynamics=Get_hash_of_initial_dynamics(Initial_Dynamics,Terms_Ordered)
print(Hash_Initial_Dynamics)

Initial_Dynamics_copy=Initial_Dynamics.copy()

t0=time_0;y0_new=y0
#................................................................

#Simplified dynamics
#................................................................
Simplified_dynamics, Simplified_times, Critical_times, Regime_results=Concatenated_simplified_dynamics(t0,time_f,y0_new,step,Initial_Dynamics_copy,Crit_Concs,Dictionary_Dynamics,Hash_Initial_Dynamics,[Switch_Predation,Switch_Burst],Terms_Ordered)

print(Critical_times)

Simplified_Solution={'time': Simplified_times, 'B_s':Simplified_dynamics[0].T , 'P_s': Simplified_dynamics[1].T }
Simplified_Solution_df=pd.DataFrame(data=Simplified_Solution)
#................................................................

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#4. ERROR OF THE SIMPLIFIED DYNAMICS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Relative Error
#--------------------------------------------------------------------
Rel_Error,Rel_Err_Bacteria,Rel_Err_Phage=Get_Relative_Error(z,Simplified_dynamics)
#--------------------------------------------------------------------

#Logratio Error
#--------------------------------------------------------------------
logratio,logratio_Bacteria,logratio_Phage=Get_logratio_Error(z,Simplified_dynamics)
#--------------------------------------------------------------------

#Total error
#--------------------------------------------------------------------
Error_bact, Error_phage=Get_Total_Error(z,Simplified_dynamics, len(z[0]))
#--------------------------------------------------------------------

#Statistics of errors
relative_bacteria_df=pd.DataFrame({'bacteria_relative':Rel_Err_Bacteria})
relative_phage_df=pd.DataFrame({'phage_relative':Rel_Err_Phage})

logratio_bacteria_df=pd.DataFrame({'bacteria_logratio':logratio_Bacteria})
logratio_phage_df=pd.DataFrame({'phage_logratio': logratio_Phage})

dataframes_errors=[relative_bacteria_df,relative_phage_df,logratio_bacteria_df,logratio_phage_df]

Statistics_Errors={}
for df_error in dataframes_errors:
    name=df_error.columns.tolist()[0]
    mean=df_error.mean().values[0]
    std=df_error.std().values[0]
    max_error=df_error.max().values[0]
    end_error=df_error.iloc[-1].values[0]
    median=df_error.median().values[0]

    #Convert dataframe to array to compute median abs deviation
    df_as_nparray=df_error.to_numpy()
    df_as_nparray=df_as_nparray.flatten()
    median_abs_deviation=scipy.stats.median_abs_deviation(df_as_nparray)
    
    Statistics_Errors[name]=[mean, std, max_error, end_error, median,median_abs_deviation]



Statistics_errors_df=pd.DataFrame(Statistics_Errors)
Statistics_errors_df['category']=['mean', 'std', 'max', 'end', 'median','mad']
Statistics_errors_df['threshold']=Active_order
Statistics_errors_df['scenario']=scenario
print(Statistics_errors_df)

Save_statistics=False
if Save_statistics==True:
	Statistics_errors_df.to_csv('../../../data/input_output_data/errors_model.csv',sep=',',mode='a',header=False, index=False)


Rel_Error_Phage_df=pd.DataFrame(Rel_Err_Phage)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#OUTPUT OF SIMULATIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
save_outputs=False

if save_outputs==True:
    tau=time_f
    weight_r=r*time_f;weight_m=m*time_f
    tau_r=1/r;tau_m=1/m
    B_crit=Critical_Bact;P_crit=Critical_Phage
    B_0=y0[0];P_0=y0[1]
    B_final=z[0][-1];P_final=z[1][-1]
    B_eq=m/(c*d);P_eq=r/d


    Outputs=[r,d,c,m,tau,weight_r,weight_m,tau_r,tau_m,B_crit,P_crit,B_0,P_0,B_final,P_final,B_eq,P_eq]
    with open('Parameters_' + scenario+ '.txt','w', encoding="utf-8") as parameters_to_save:
        for item in Outputs:
            parameters_to_save.write( str(item) + "\n")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#Plots
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#3.1. Configuration of plots
#====================================================================
#Path to save figure
Output_Path='../../../results/Plots/plots_paper/'

#Fontsizes
size_axis=7;size_ticks=6;size_title=5

#Figure Size
cm = 1/2.54  # centimeters in inches
Width=8*cm;Height=4*cm #Width and height of plots

#Gridspec parameters
Rows=1;Cols=2

#Colors
cmap='RdBu';cmap_pieces= plt.get_cmap(cmap)
Phage_color=cmap_pieces(0.1);Bacteria_color=cmap_pieces(0.9)

#Extensions to save
Extensions=['.svg','.pdf']

#Linewidth
width_line=1
#====================================================================


#Plot relative error
#====================================================================
fig=figure(figsize=(Width, Height), dpi=300)


#Select type of error
#--------------------------------------------------
type_error='relative' #log-ratio, relative

if type_error=='log-ratio':
    error_bacteria=logratio_Bacteria
    error_phage=logratio_Phage
    error=logratio

    max_error=np.max(error); time_max_error=time[error.index(max_error)]
    end_error=error[-1]
    mean_error=np.mean(error)
    
    ylabel_error='Logratio error'
    Name_Figure='Logratio_Error_'+str(scenario)


#    ylimits= [ [-0.5,70], [-0.1,5] , [-1,90], [-0.5,70] ]

    ylimits= [ 0.90, 0.1 , 0.90, 0.90]
    
elif type_error=='relative':
    error_bacteria=Rel_Err_Bacteria
    error_phage=Rel_Err_Phage   
    error=Rel_Error

    max_error=np.max(error); time_max_error=time[error.index(max_error)]
    end_error=error[-1]
    mean_error=np.mean(error)

    
    ylabel_error='Relative error'
    Name_Figure='Relative_Error_'+str(scenario)

#    ylimits= [ [-0.5,55], [-0.1,5] , [-1,90], [-0.5,70] ]

    ylimits= [ 0.9, 0.05 , 1.5, 0.5 ]
#--------------------------------------------------

#Percentage error representation
#--------------------------------------------------
percentage=False
if percentage==True:
    error_bacteria=[100*i for i in error_bacteria]
    error_phage=[100*j for j in error_phage]
    error=[100*k for k in error]
    max_error=100*max_error

    ylimits=[ylim*100 for ylim in  ylimits]
    ylabel_error= ylabel_error  + ' (%)'
#--------------------------------------------------

try:
    plt.plot(time,error_bacteria,color=Bacteria_color,linewidth=1,label='Prey')
    plt.plot(time,error_phage,color=Phage_color,linewidth=1,label='Predator')
    plt.plot(time,error,color='k',linewidth=1,label='Mean')
    
except ValueError:
    plt.plot(Simplified_times,error_bacteria,color=Bacteria_color,linewidth=1,label='Prey')
    plt.plot(Simplified_times,error_phage,color=Phage_color,linewidth=1,label='Predator')
    plt.plot(Simplified_times,error,color='k',linewidth=1,label='Mean')

for crit_time in Critical_times:
    plt.axvline(x=crit_time, linestyle='dotted', color='k', linewidth=1)

    
plt.xlabel('Time (h)',fontsize=size_axis);
plt.ylabel(ylabel_error,fontsize=size_axis)

if scenario=='growth':
    plt.ylim(top=ylimits[0])
elif scenario=='decay':
    plt.ylim(top=ylimits[1])
elif scenario=='growth-decay_regimes':
    plt.ylim(top=ylimits[2])
elif scenario=='no-growth-no-decay':
    plt.ylim(top=ylimits[3])
    
if time_f<=20:
    step_ticks=2
elif time_f>20 and time_f<=100:
    step_ticks=10
elif time_f>100 and time_f<=500:
    step_ticks=50
elif time_f>500 and time_f<=1000:
    step_ticks=100
elif time_f>1000:
    step_ticks=250
else:
    step_ticks=int(0.1*time_f)

xticks_labels=[label for label in range(0,time_f+1,step_ticks)]

plt.xticks(xticks_labels, fontsize=size_ticks)

#y logscale
#-------------------
# plt.yscale('log');
# plt.ylim([5e-9,30])
# yticks_labels=[1e-8, 1e-5, 1e-2, 10]
# plt.yticks(yticks_labels, fontsize=size_ticks)
#-------------------

plt.yticks(fontsize=size_ticks)

if scenario=='growth':
    plt.legend(loc='best',fontsize=size_ticks,frameon=False)

for ext in Extensions:
    plt.savefig(Output_Path+Name_Figure+ext,dpi=300)

plt.show()
#====================================================================

'''
#Plot phase diagram
#====================================================================
#Figure and gridspec
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
fig=figure(figsize=(Width, Height), dpi=300)


gs=gridspec.GridSpec(1,2)
gs.update(left=0.1,right=0.95,bottom=0.08,top=0.97,wspace=0.1,hspace=0.0)

ax_00=plt.subplot(gs[0,0])
#-----------------------------------------------------

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Plot Figure
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
plt.plot(z[1].T,z[0].T,color='k',linewidth=1,label='Full')
plt.scatter(z[1][0], z[0][0],c='gray')
plt.scatter(z[1][-1], z[0][-1],c='gray')
matplotlib.pyplot.text(z[1][0], z[0][0], "$t_0$",fontsize=size_ticks)
matplotlib.pyplot.text(z[1][-1], z[0][-1],"$t_f$",fontsize=size_ticks)

Max_value_phage=np.max(Simplified_dynamics[1].T)
Min_value_phage=np.min(Simplified_dynamics[1].T)

Max_value_bact=np.max(Simplified_dynamics[0].T)
Min_value_bact=np.min(Simplified_dynamics[0].T)

plt.hlines(Crit_Concs[1],xmin=0.5*Min_value_phage,xmax=15*Max_value_phage,color=Bacteria_color,linewidth=1,linestyle='-')
plt.vlines(Crit_Concs[0],ymin=0.5*Min_value_bact,ymax=20*Max_value_bact,color=Phage_color,linewidth=1,linestyle='-')

plt.plot(Simplified_dynamics[1], Simplified_dynamics[0],color='gray',linestyle='dashed',dashes=(5,5),linewidth=1,label='Simplified')
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Labels, ticks, and legend
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
plt.xlabel('Predator concentration $\mathrm{(ml^{-1})}$',fontsize=size_axis)
plt.ylabel('Prey concentration $\mathrm{(ml^{-1})}$',fontsize=size_axis)

plt.xscale('log');plt.yscale('log')

plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks) 
plt.minorticks_off()

plt.legend(loc='best',fontsize=size_ticks,frameon=False)

xlimits=plt.gca().get_xlim()    
ylimits=plt.gca().get_ylim()

if scenario=='growth-decay':
    plt.xlim([8852.684347833512,552195289.6473755])
    plt.ylim([58.17505512731038,4979516.565870757])

    plt.xlim([3956.8362369579113, 653621046.5628147])
    plt.ylim([26.002188481729828, 5894139.225268255])



sns.despine(top=True, right=True, left=False, bottom=False)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Save figure
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
Name_Figure='Phase_Diagrams_'+str(scenario)

for ext in Extensions:
    plt.savefig(Output_Path+Name_Figure+ext,dpi=300)

plt.show()
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

if scenario=='growth':
    fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True, gridspec_kw={'height_ratios':[3,1]}, figsize=(Width,Height), dpi=300)
    
    fig.subplots_adjust(hspace=0.2)  # adjust space between axes

    ax1.plot(z[1].T,z[0].T,color='k',linewidth=1,label='Full')
    ax1.plot(Simplified_dynamics[1], Simplified_dynamics[0],color='gray',linestyle='dashed',dashes=(5,5),linewidth=1,label='Simplified')
    ax1.scatter(z[1][0], z[0][0],c='gray')
    
    ax1.text(z[1][0], z[0][0], "$t_0$",fontsize=size_ticks)
    ax1.set_xscale('log');ax1.set_yscale('log')
    ax1.set_ylim(500, 1e7)
    ax1.vlines(Crit_Concs[0],ymin=0.5*Min_value_bact,ymax=20*Max_value_bact,color=Phage_color,linewidth=1,linestyle='-')
    ax1.hlines(Crit_Concs[1],xmin=0.5*Min_value_phage,xmax=15*Max_value_phage,color=Bacteria_color,linewidth=1,linestyle='-')

    ax1.set_yticks([1e3, 1e6]) 
    ax1.tick_params(axis='both', which='major', labelsize=size_ticks)
    ax1.set_ylabel('Prey concentration, B $\mathrm{(ml^{-1})}$',fontsize=size_axis)
    plt.minorticks_off()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.xaxis.set_ticks_position('none')
    ax1.minorticks_off()
    
    ax2.plot(z[1].T,z[0].T,color='k',linewidth=1,label='Full')
    ax2.plot(Simplified_dynamics[1], Simplified_dynamics[0],color='gray',linestyle='dashed',dashes=(5,5),linewidth=1,label='Simplified')
    
    ax2.scatter(z[1][-1], z[0][-1],c='gray')
    matplotlib.pyplot.text(z[1][-1], z[0][-1],"$t_f$",fontsize=size_ticks)
    ax2.text(z[1][-1], z[0][-1],"$t_f$",fontsize=size_ticks)    

    ax2.set_xscale('log');ax2.set_yscale('log')
    
    ax2.set_ylim(4.0019028108328327e-16, 1e-8)
    ax2.vlines(Crit_Concs[0],ymin=0.5*Min_value_bact,ymax=20*Max_value_bact,color=Phage_color,linewidth=1,linestyle='-')
    ax2.hlines(Crit_Concs[1],xmin=0.5*Min_value_phage,xmax=15*Max_value_phage,color=Bacteria_color,linewidth=1,linestyle='-')

    ax2.set_yticks([1e-10])
    
    ax2.set_xticks([1e5, 1e7,1e9]) 
    ax2.tick_params(axis='both', which='major', labelsize=size_ticks)
    ax2.set_xlabel('Predator concentration, P $\mathrm{(ml^{-1})}$',fontsize=size_axis)
    plt.minorticks_off()


    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.legend(loc='best',fontsize=size_ticks,frameon=False)
    
#    ax1.tick_params(labeltop=False)# don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    
    d = .3  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)

    ax1.plot([0], [0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0], [1], transform=ax2.transAxes, **kwargs)


    Name_Figure='Broken_Phase_Diagram_'+str(scenario)

    for ext in Extensions:
        plt.savefig(Output_Path+Name_Figure+ext,dpi=300)

    plt.show()
#====================================================================

#Plot Full Dynamics
#====================================================================

#Figure and gridspec
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
fig=figure(figsize=(Width, Height), dpi=300)

Gridspec_layout=False
if Gridspec_layout==True:
    gs=gridspec.GridSpec(Rows,Cols,width_ratios=[1,0.04])
    gs.update(left=0.1,right=0.7,bottom=0.08,top=0.97,wspace=0.1,hspace=0.1)
    ax_00=plt.subplot(gs[0,0])
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
plt.plot(time,z[0].T,color=Bacteria_color,linewidth=1,label='Prey')
plt.plot(time,z[1].T,color=Phage_color,linewidth=1,label='Predator')

plt.plot(Simplified_times,Simplified_dynamics[0].T,color='k',linestyle='dashed',dashes=(5,5),linewidth=1,label='Simplified')

plt.plot(Simplified_times,Simplified_dynamics[1].T,color='k',linestyle='dashed',dashes=(5,5),linewidth=1)

for crit_time in Critical_times:
    plt.axvline(x=crit_time, linestyle='dotted', color='k', linewidth=1)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Labels, axes, ticks, and legend
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
plt.ylabel('Concentration $\mathrm{(ml^{-1})}$',fontsize=size_axis)
#plt.ylabel('$A_1$',fontsize=size_axis)
#plt.xlabel('Time (h)',fontsize=size_axis)
plt.xlabel('Time',fontsize=size_axis)

#x,y limits
plt.xlim([time[0],time[-1]])
plt.yscale('log')
plt.ylim([1, 10*np.max(z)])

time_0=int(time_0);time_f=int(time_f)

if time_f<=20:
    step_ticks=2
elif time_f>20 and time_f<=100:
    step_ticks=10
elif time_f>100 and time_f<=500:
    step_ticks=50
elif time_f>500 and time_f<=1000:
    step_ticks=100
elif time_f>1000:
    step_ticks=250
else:
    step_ticks=int(0.1*time_f)

xticks_labels=[label for label in range(0,time_f+1,step_ticks)]
     
plt.xticks(xticks_labels)
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks) 

plt.legend(loc='best',fontsize=size_ticks,frameon=False)
sns.despine(top=True, right=True, left=False, bottom=False)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Save figure
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
Name_Figure='Dynamics_tau_'+str(time_f)

for ext in Extensions:
    plt.savefig(Output_Path+Name_Figure+ext,dpi=300)
plt.show()
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
#====================================================================


#Plot sum of mechanisms
#====================================================================

#Figure and gridspec
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
fig=figure(figsize=(Width, Height), dpi=300)

Gridspec_layout=False
if Gridspec_layout==True:
    gs=gridspec.GridSpec(Rows,Cols,width_ratios=[1,0.04])
    gs.update(left=0.1,right=0.7,bottom=0.08,top=0.97,wspace=0.1,hspace=0.1)
    ax_00=plt.subplot(gs[0,0])
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Plot
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
plt.plot(time,Total_Active_Bioprocesses,color='k',linewidth=width_line,linestyle='-')

for crit_time in Critical_times:
    plt.axvline(x=crit_time, linestyle='dotted',color='k', linewidth=1)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Labels, axes, ticks, and legend
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
if time_f<=20:
    step_ticks=2
elif time_f>20 and time_f<=100:
    step_ticks=10
elif time_f>100 and time_f<=500:
    step_ticks=50
elif time_f>500 and time_f<=1000:
    step_ticks=100
elif time_f>1000:
    step_ticks=250
else:
    step_ticks=int(0.1*time_f)

xticks_labels=[label for label in range(0,time_f+1,step_ticks)]

plt.ylabel('Number of active terms',fontsize=size_axis)
plt.xlabel('Time (h)',fontsize=size_axis)

plt.xlim([time_0,time_f]);plt.ylim([0,4.5])

plt.xticks(xticks_labels)
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)

plt.legend(loc='best',fontsize=size_ticks,frameon=False)

sns.despine(top=True, right=True, left=False, bottom=False)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Save figure
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
Name_Fig='Sum_of_Mechs'+str(scenario)

[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]
plt.show()
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
#====================================================================

#Figure weights
#====================================================================
fig=figure(figsize=(Width, Height), dpi=300)
Rows_weights=4
Cols_weights=1

growth_weights=Bioprocesses_df['Growth'].tolist()
predation_weights=Bioprocesses_df['Predation'].tolist()
burst_weights=Bioprocesses_df['Burst'].tolist()
decay_weights=Bioprocesses_df['Decay'].tolist()


gs=gridspec.GridSpec(Rows_weights,Cols_weights)
gs.update(left=0.1,right=0.95,bottom=0.08,top=0.97,wspace=0.1,hspace=0.0)


Time_Bioprocesses_df=Bioprocesses_df.index.values.tolist()
Bioprocesses_df.set_index(time,inplace=True)

cols_0 = list(Bioprocesses_df.columns.values)
Bioprocesses_df = Bioprocesses_df[['Burst', 'Decay', 'Growth', 'Predation']]

Bioprocesses_df_transp = Bioprocesses_df.transpose()
Max_val_0=Bioprocesses_df.max();Max_val=np.max(Bioprocesses_df.max())
Min_val_0=Bioprocesses_df.min();Min_val=np.min(Bioprocesses_df.min())
y_lim_max=math.ceil(math.log(Max_val,10))

Constant_Processes=Bioprocesses_df[['Decay','Growth']]
Min_Constant_Processes=np.min(Constant_Processes.min())

y_lim_min=math.floor(math.log(Min_Constant_Processes,10))
if y_lim_min>0:
    y_lim_min=0

xticks_pos=[position for position in np.arange(Time_Bioprocesses_df[0],Time_Bioprocesses_df[-1]+1,step_ticks*Resolution)]

order_step=2
yticks_pos=[10**(y_lim_min),10**(math.floor(0.5*(y_lim_min+y_lim_max))),10**(y_lim_max)] 

ax_00=plt.subplot(gs[0,0])
#-----------------------------------------------------
plt.plot(time,growth_weights,color=Bacteria_color)
plt.axhline(y=Active_order, linestyle='dotted', color='gray', linewidth=1)

plt.ylabel('w$_{g}$',fontsize=size_axis)
plt.yscale('log')

plt.ylim([(10**(y_lim_min-1)), 9*(10**y_lim_max)])

plt.minorticks_off()
ax_00.set_yticks(yticks_pos)
yticks = ax_00.yaxis.get_major_ticks()
plt.yticks(fontsize=size_ticks)
yticks[1].label1.set_visible(False)

plt.xticks(xticks_labels)
plt.xticks(fontsize=size_ticks)
#-----------------------------------------------------

ax_10=plt.subplot(gs[1,0],sharex=ax_00)
#-----------------------------------------------------
plt.plot(time,predation_weights,color=Bacteria_color)
plt.axhline(y=1, linestyle='dotted', color='gray', linewidth=1)

plt.ylabel('w$_{p}$',fontsize=size_axis)
plt.yscale('log')
plt.ylim([(10**(y_lim_min-1)), 9*(10**y_lim_max)])

plt.minorticks_off()
ax_10.set_yticks(yticks_pos)
yticks = ax_10.yaxis.get_major_ticks()
plt.yticks(fontsize=size_ticks)
yticks[1].label1.set_visible(False)

plt.xticks(xticks_labels)
plt.xticks(fontsize=size_ticks)
#-----------------------------------------------------

ax_20=plt.subplot(gs[2,0],sharex=ax_00)
#-----------------------------------------------------
plt.plot(time,burst_weights,color=Phage_color)
plt.axhline(y=1, linestyle='dotted', color='gray', linewidth=1)

plt.ylabel('w$_{b}$',fontsize=size_axis)
plt.yscale('log')
plt.ylim([(10**(y_lim_min-1)), 9*(10**y_lim_max)])

plt.minorticks_off()
ax_20.set_yticks(yticks_pos)
yticks = ax_20.yaxis.get_major_ticks()
plt.yticks(fontsize=size_ticks)
yticks[1].label1.set_visible(False)

plt.xticks(xticks_labels)
plt.xticks(fontsize=size_ticks)
#-----------------------------------------------------

ax_30=plt.subplot(gs[3,0],sharex=ax_00)
#-----------------------------------------------------
plt.plot(time,decay_weights,color=Phage_color)
plt.axhline(y=1, linestyle='dotted', color='gray', linewidth=1)

plt.ylabel('w$_{d}$',fontsize=size_axis)
plt.yscale('log')
plt.ylim([(10**(y_lim_min-1)), 9*(10**y_lim_max)])

plt.minorticks_off()
yticks = ax_30.yaxis.get_major_ticks()
ax_30.set_yticks(yticks_pos)
plt.yticks(fontsize=size_ticks)
yticks[1].label1.set_visible(False)

plt.xticks(xticks_labels)
plt.xticks(fontsize=size_ticks)
#-----------------------------------------------------

plt.xticks(fontsize=size_ticks);
plt.xlabel('Time (h)',fontsize=size_axis)


Name_Fig='Weights_'+str(scenario)
[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]
            
plt.show()
#====================================================================

#Heatmap
#====================================================================

#Figure and gridspec
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
fig=figure(figsize=(Width, Height), dpi=300)

Gridspec_layout=False
if Gridspec_layout==True:
    gs=gridspec.GridSpec(Rows,Cols,width_ratios=[1,0.04])
    gs.update(left=0.1,right=0.7,bottom=0.08,top=0.97,wspace=0.1,hspace=0.1)
    ax_00=plt.subplot(gs[0,0])
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::


#3.4.1. Ticks, labels, axes, and legend
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if Dominant_Timescale=='r':

if time_f<=20:

    ax_heatmap=sns.heatmap(Bioprocesses_df_transp,vmin=0.0,vmax=Max_val,cmap="Greys",cbar=True)
    
    step_ticks=2 
    xticks_pos_heatmap=[position for position in np.arange(Time_Bioprocesses_df[0],Time_Bioprocesses_df[-1]+1,step_ticks*Resolution)]
    xticks_labels_heatmap=[label for label in np.arange(time_0,time_f+1,step_ticks)]

elif time_f>20 and time_f<=100:
    
    ax_heatmap=sns.heatmap(Bioprocesses_df_transp,vmin=0.0,vmax=Max_val,cmap="Greys",cbar=True)
    
    step_ticks=10
    xticks_pos_heatmap=[position for position in np.arange(Time_Bioprocesses_df[0],Time_Bioprocesses_df[-1]+1,step_ticks*Resolution)]
    xticks_labels_heatmap=[label for label in np.arange(time_0,time_f+1,step_ticks)]
    
elif time_f>100 and time_f<=500:

    #Sample the dataframe
    #--------------------------------------------------------------
    Light_Heatmap_Vector=time[0::100]
    Light_Bioprocesses_df=Bioprocesses_df.loc[Light_Heatmap_Vector]
    Light_Bioprocesses_df_transp=Light_Bioprocesses_df.transpose()
    Light_time=Light_Bioprocesses_df.index.values.tolist()
    Light_time=np.asarray(Light_time)
    #--------------------------------------------------------------

    step_ticks=50
    xticks_pos_heatmap=[position for position in np.arange(Light_time[0],Light_time[-1]+1,step_ticks*Resolution/100)]

    xticks_labels_heatmap=[int(label) for label in np.arange(time_0,time_f+1,step_ticks)]

    ax_heatmap=sns.heatmap(Light_Bioprocesses_df_transp,vmin=Min_val,vmax=Max_val,cmap="Greys",cbar=False)
    
elif time_f>500:

    #Sample the dataframe
    #--------------------------------------------------------------
    Light_Heatmap_Vector=time[0::100]
    Light_Bioprocesses_df=Bioprocesses_df.loc[Light_Heatmap_Vector]
    Light_Bioprocesses_df_transp=Light_Bioprocesses_df.transpose()
    Light_time=Light_Bioprocesses_df.index.values.tolist()
    Light_time=np.asarray(Light_time)
    #--------------------------------------------------------------
    
    step_ticks=100

    xticks_pos_heatmap=[position for position in np.arange(Light_time[0],Light_time[-1]+1,step_ticks*Resolution/100)]

    xticks_labels_heatmap=[int(label) for label in np.arange(time_0,time_f+1,step_ticks)]

    ax_heatmap=sns.heatmap(Light_Bioprocesses_df_transp,vmin=Min_val,vmax=Max_val,cmap="Greys",cbar=False)


ax_heatmap.set_xticks(xticks_pos_heatmap)
ax_heatmap.set_xticklabels(xticks_labels_heatmap,fontsize=size_ticks,rotation='horizontal')


plt.yticks(fontsize=size_ticks)
ax_heatmap.tick_params(axis=u'y', which=u'both',length=0)
plt.xlabel('Time (h)',fontsize=size_axis)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Save Figure
#:::::::::::::::::::::::::::::::::::::::::::::::::::::
Name_Fig='Terms_'+str(scenario)
[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]
            
plt.show()
#:::::::::::::::::::::::::::::::::::::::::::::::::::::
#====================================================================

'''
