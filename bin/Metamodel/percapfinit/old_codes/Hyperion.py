#2021-09-07


#Import libraries
#+++++++++++++++++++++++++++++++++++++++++++++++++++
from Metamodel_Functions import *
import numpy as np
from random import shuffle
import matplotlib
import matplotlib.gridspec as gridspec
import time
import matplotlib.pyplot as plt
import seaborn as sns
#+++++++++++++++++++++++++++++++++++++++++++++++++++


#1. Latin Hypercube Sampling
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Matt's parameters
Delta_Interval=0.1

#Growth rates
mean_r=0.9
range_r=[(1-Delta_Interval)*mean_r,(1+Delta_Interval)*mean_r]

#Infection rates
mean_d=3e-8
range_d=[(1-Delta_Interval)*mean_d,(1+Delta_Interval)*mean_d]

#Burst sizes
mean_c=150
range_c=[(1-Delta_Interval)*mean_c,(1+Delta_Interval)*mean_c]

#Phage decay rates
mean_m=2.8e-3
range_m=[(1-Delta_Interval)*mean_m,(1+Delta_Interval)*mean_m]

#Carrying capacities
mean_K=1e11
range_K=[(1-Delta_Interval)*mean_K,(1+Delta_Interval)*mean_K]

#Sampling Points
N=10

Ranges_Parameters={'r':range_r,'d':range_d,'c':range_c,'m':range_m,'K':range_K}
Samples_Parameters=LHS(Ranges_Parameters, N, 3333)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. Solve equations
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Initial Conditions
#================
B0=1000
T0=10000
Initial_Conditions=[B0,T0]
#================

#Time and resolution
#======================================
time_0=0;time_f=14
Resolution=100
Steps=time_f*Resolution+1
time=np.linspace(time_0,time_f,Steps)
step=0.01
#====================================== 

#Run over samplings of LHS and solve L-V equations
#============================================================
Solutions={};Thresholds_Dict={};Threshold_Times={}

for sampling in Samples_Parameters:

    #Extract parameters from LHS
    #.......................................................................
    Parameters=Samples_Parameters[sampling]

    r=Parameters['r'];K=Parameters['K'];m=Parameters['m'];
    c=Parameters['c'];d=Parameters['d']
    #.......................................................................

    #Solve and extract concentrations
    #.......................................................................
    Solution=Solve_Lotka_Volterra(Parameters,Initial_Conditions, time, step)
    Solutions[sampling]=Solution
    z=Solutions[sampling].sol(time)
    Concentrations={'Bacteria':z[0], 'Phage':z[1]}
    #.......................................................................


    #Extract times at which thresholds are crossed. Save times+thresholds in dicts
    #............................................................................
    Thresholds_Sampling={'Bacteria':{'B_growth': K,'B_burst': r/(c*d)},\
                                   'Phage':{'T_infection':(r/d)} }
    for Thing in Thresholds_Sampling:
        for Term in Thresholds_Sampling[Thing]:

            Threshold=Thresholds_Sampling[Thing][Term]
            Thing_Concentration=Concentrations[Thing]
            print("hola")
            print(type(Threshold))
            print(type(Thing_Concentration))
            Threshold_crossings=np.where(np.sign(Thing_Concentration[:-1] - Threshold) != np.sign(Thing_Concentration[1:]-Threshold) )[0] + 1

            #Save thresholds in dict
            #---------------------------------------------------------------
            try:
                Thresholds_Dict[Term][sampling]=Thresholds_Sampling[Thing][Term]
            except KeyError:
                Thresholds_Dict[Term]={sampling:Thresholds_Sampling[Thing][Term]}
            #---------------------------------------------------------------

            #Save times of thresholds in dict
            #---------------------------------------------------------------
            if len(Threshold_crossings>0):#skip not found critical times
                Threshold_Times_Term=[time[indx] for indx in Threshold_crossings]
                
                try:
                    Threshold_Times[Term][sampling]=Threshold_Times_Term
                    
                except KeyError:
                    Threshold_Times[Term]=\
                            {sampling:Threshold_Times_Term}
                else:
                    continue
            #---------------------------------------------------------------
#===============================================================================
print(Threshold_Times)

#Mean concentrations and times of thresholds
#===============================================================================

#Mean times
#-----------------------------------------------------------
Mean_Threshold_Times={};Std_Threshold_Times={}
for Term in Threshold_Times:
    
    Values_Term=Threshold_Times[Term].values()
    Values_Term=list(Values_Term)
    print(Term)
    print(Values_Term)
    Mean_Values=[];Std_Values=[]
    print(len(Values_Term))
    print(len(Values_Term[0]))
    
#    for i in range(len(Values_Term)):
    for i in range(len(Values_Term[0])):
        print(i)
        l = [item[i] for item in Values_Term]
        print(l)
        Mean_Values.append(np.mean(l))
        Std_Values.append(np.std(l))

    Mean_Threshold_Times[Term]=Mean_Values
    Std_Threshold_Times[Term]=Std_Values
#-----------------------------------------------------------

#Mean Concentrations and STDs over LHS
#------------------------------------------------------------------------------

#Convert dicts to dataframes
#........................................................................
Thresholds_df=pd.DataFrame(Thresholds_Dict)

Bacterial_Concentrations={sampling: Solutions[sampling].sol(time)[0] for sampling in Solutions}
Phage_Concentrations={sampling: Solutions[sampling].sol(time)[1] for sampling in Solutions}

Bacterial_Concentrations_df=pd.DataFrame(Bacterial_Concentrations)
Phage_Concentrations_df=pd.DataFrame(Phage_Concentrations)
#........................................................................

#Means and stds
#........................................................................
Mean_Bacterial_Concentrations=Bacterial_Concentrations_df.mean(axis=1)
Std_Bacterial_Concentrations=Bacterial_Concentrations_df.std(axis=1)

Bacterial_Under_line=(Mean_Bacterial_Concentrations - 0.5*Std_Bacterial_Concentrations)
Bacterial_Over_line =(Mean_Bacterial_Concentrations + 0.5*Std_Bacterial_Concentrations)

Mean_Phage_Concentrations=Phage_Concentrations_df.mean(axis=1)
Std_Phage_Concentrations=Phage_Concentrations_df.std(axis=1)

Phage_Under_line=(Mean_Phage_Concentrations - 0.5*Std_Phage_Concentrations)
Phage_Over_line =(Mean_Phage_Concentrations + 0.5*Std_Phage_Concentrations)
#........................................................................
#===============================================================================



#Solution for a single sampling
#===========================================================================
#Dictionary to save the values of individual terms
Bioprocesses={'Growth':[], 'Infection':[], 'Burst':[], 'Decay':[]}

Colors={'B_growth':'g', 'B_burst':'r', 'T_infection':'b'}

Sampling_Point=0
Parameters=Samples_Parameters[Sampling_Point]

r=Parameters['r'];K=Parameters['K'];m=Parameters['m'];
c=Parameters['c'];d=Parameters['d']

#Critical concentrations
# Thresholds_Dict={'Bacteria':{'B_growth':{K}, 'B_burst':{m/(c*d)}},
#                                       'Phage':{'T_infection':{r/d} }}

Thresholds_Dict={'Bacteria':{'B_growth':{K}, 'B_burst':{r/(c*d)}},
                                      'Phage':{'T_infection':{r/d} }}

z=Solutions[Sampling_Point].sol(time)

Concentrations={'Bacteria':z[0], 'Phage':z[1]}


Threshold_Times={}
for Thing in Thresholds_Dict:
    for Term in Thresholds_Dict[Thing]:
        Threshold=list(Thresholds_Dict[Thing][Term])[0]
        Concentration=Concentrations[Thing]
        Threshold_crossings=np.where(np.sign(Concentration[:-1] - Threshold) != np.sign(Concentration[1:]-Threshold) )[0] + 1

        if len(Threshold_crossings>0):
            Threshold_Times[Term]=Threshold_crossings
        else:
            continue


for (Bacteria,Phage) in zip(z[0],z[1]):

    Growth_t=(r*Bacteria)
    Normalized_Growth_t=Growth_t/(Bacteria*r)
    
    Infection_t=(d*Bacteria*Phage)
    Normalized_Infection_t=Infection_t/(Bacteria*r)
    
    Burst_t=(c*d*Bacteria*Phage)
    Normalized_Burst_t=Burst_t/(Phage*r)
    
    Decay_t=(m*Phage)
    Normalized_Decay_t=Decay_t/(Phage*r)


    #Append values of individual terms
    #.........................................................
    Bioprocesses['Growth'].append(Normalized_Growth_t)
    Bioprocesses['Infection'].append(Normalized_Infection_t)
    Bioprocesses['Burst'].append(Normalized_Burst_t)
    Bioprocesses['Decay'].append(Normalized_Decay_t) 
    #.........................................................
#===============================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Path to save
#================================
path_plot='../../../results/Plots/'
#================================

#1.4. Equilibrium concentrations
#==============================
B_eq=m/(c*d)
T_eq=(r/d)*(1-B_eq/K)
#============================== 

#3.1. Convert Information to Dataframes
#===========================================
Dataframe=pd.DataFrame(data=Bioprocesses)


Relative_Tension_Growth=Dataframe['Growth']/Dataframe['Burst']
Relative_Tension_Infect=Dataframe['Infection']/Dataframe['Decay']
Diff_Tension_Growth=Dataframe['Growth']-Dataframe['Burst']
Diff_Tension_Infect=Dataframe['Infection']-Dataframe['Decay']
Bacteria_Relative_Equilibrium=z[0]/B_eq
Phage_Relative_Equilibrium=z[1]/T_eq
Dataframe_transp = Dataframe.transpose()

Max_val_0=Dataframe.max()
Max_val=np.max(Dataframe.max())
Min_val_0=Dataframe.min()
Min_val=np.min(Dataframe.min())
#===========================================

#Parameters and fontsizes
#========================
size_axis=18
size_ticks=16
size_letter=15
Letters=['a','b']
#=======================


#First plot
#=================================================================================
# #Gridspec
# #----------------------------------------------------------------------------
# fig=plt.figure(figsize=(11, 15))
# gs=gridspec.GridSpec(1, 1, height_ratios=[1])
# gs.update(left=0.07,right=0.91,bottom=0.08,top=0.97,wspace=0.2,hspace=0.2)
# #----------------------------------------------------------------------------

# #Concentrations
# #----------------------------------------------------------------------------
# ax_00=plt.subplot(gs[0,0])


# plt.plot(time,Mean_Bacterial_Concentrations, linewidth=2) #mean curve.
# plt.fill_between(Std_Bacterial_Concentrations.index, Bacterial_Under_line, Bacterial_Over_line, color='b', alpha=.1)

# plt.plot(time,Mean_Phage_Concentrations, linewidth=2,color='r') #mean curve.
# plt.fill_between(Std_Phage_Concentrations.index, Phage_Under_line, Phage_Over_line, color='r', alpha=.1)

# plt.axvline(x=Mean_Threshold_Times['B_burst'],color='k',linewidth=3,linestyle='--')

# plt.axvspan(Mean_Threshold_Times['B_burst'][0]-0.5*Std_Threshold_Times['B_burst'][0], Mean_Threshold_Times['B_burst']+0.5*Std_Threshold_Times['B_burst'][0],alpha=0.5, color='blue')

# Smallest_Time=Mean_Threshold_Times['T_infection'][0]-0.5*Std_Threshold_Times['T_infection'][0]
# Smallest_Measure_Time=Smallest_Time-1

# s = np.random.uniform(Smallest_Measure_Time,Smallest_Time,10)
# print(s)
# print(time)
# plt.axvspan(Smallest_Measure_Time, Smallest_Time,alpha=0.5, color='gray')
# for element in s:
#     plt.axvline(x=element,color='k',linestyle='--',linewidth=0.5)

# # for value in s:
# #     print value
    
# #     Closest_Time=min(range(len(time)), key=lambda i: abs(time[i]-value))

# #     print(Closest_Time)
# #     print(time[Closest_Time])
# #     print(s[Closest_Time])
# #     print(time[6])
# #     print('\n')

# # plt.axvline(x=Mean_Threshold_Times['T_infection'],color='k',linewidth=3)
# # plt.axvspan(Mean_Threshold_Times['T_infection'][0]-0.5*Std_Threshold_Times['T_infection'][0], Mean_Threshold_Times['T_infection']+0.5*Std_Threshold_Times['T_infection'][0],alpha=0.5, color='red')

# plt.ylabel(r'Concentration $(ml^{-1})$',fontsize=size_axis)
# plt.xlabel(r'Time $(h)$',fontsize=size_axis)

# ax_00.set_ylim([1,1e10])
# ax_00.set_yscale('log')
#----------------------------------------------------------------------------

#Boxplots
#----------------------------------------------------------------------------
#ax_00=plt.subplot(gs[1,0])

#boxplot = Thresholds_df.boxplot(column=['B_burst', 'T_infection'])

#plt.ylabel(r'Concentration $(ml^{-1})$',fontsize=size_axis)
#plt.xlabel(r'Time $(h)$',fontsize=size_axis)

#plt.yscale('log')
#----------------------------------------------------------------------------

#Save and show
#----------------------------------------------------------------------------
# plt.savefig(path_plot + 'Hyperion_LHS_'+ '.pdf',dpi=300)
# plt.savefig(path_plot + 'Hyperion_LHS_'+ '.png',dpi=300)
# plt.show()
#----------------------------------------------------------------------------

#=================================================================================




#3.2. Configure gridspec plot
#===============================================
#Plots with tensions
fig=plt.figure(figsize=(40, 30))
#gs=gridspec.GridSpec(4, 2, width_ratios=[1,0.04])
gs=gridspec.GridSpec(2, 2, width_ratios=[1,0.04])
gs.update(left=0.07,right=0.85,bottom=0.08,top=0.97,wspace=0.12,hspace=0.02)
#===============================================

#3.3. (0,0) - Lotka-Volterra Solutions
#===============================================
ax_00=plt.subplot(gs[0,0])
plt.plot(time, z.T,linewidth=3)



#Ticks, legend, axes
#::::::::::::::::::::::::::::::::::::::::::::::: 
plt.xlim([time_0, time_f])
ax_00.set_yscale('log')
plt.yticks(fontsize=size_ticks)
plt.ylabel(r'Concentration $(ml^{-1})$',fontsize=size_axis)
plt.legend([r'E. Coli_Eq', 'T4_Eq', r'E. Coli', 'T4' ],loc='upper right')
ax_00.set_ylim([1,1e10])

print(Threshold_Times)
for Thing in Threshold_Times:
     for time_index in Threshold_Times[Thing]:
        plt.axvline(x=time[time_index],color=Colors[Thing],linewidth=3)
        
ax_00.axes.xaxis.set_ticks([])


#Letter for caption
#---------------------------------------
Lim_x_up=plt.gca().get_xlim()
Lim_y_up=plt.gca().get_ylim()
#---------------------------------------
#:::::::::::::::::::::::::::::::::::::::::::::::
#===============================================  

#3.4. (1,0) - Heatmap 
#========================================================================
ax_01 =plt.subplot(gs[1,0])
sns.heatmap(Dataframe_transp,vmin=Min_val,vmax=Max_val,cmap="Blues",cbar=False)

#3.4.1. Ticks, labels, and axes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Time_Dataframe=Dataframe.index.values.tolist()


Tickstep=1                                                                 
xticks_pos=[position for position in range(Time_Dataframe[0],Time_Dataframe[-1]+1,Tickstep*Resolution)]
xticks_labels=[label for label in range(time_0,time_f+1,Tickstep)]
ax_01.set_xticks(xticks_pos)
ax_01.set_xticklabels(xticks_labels)   
plt.xticks(rotation='horizontal',fontsize=size_ticks)

plt.yticks(fontsize=size_ticks)


plt.xlabel(r'time $(h)$',fontsize=size_axis)
Lim_x_up=plt.gca().get_xlim()
Lim_y_up=plt.gca().get_ylim()            
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#========================================================================

#3.5. (1,1) - Colorbar
#===============================================
# colorbar=plt.subplot(gs[1,1])
# colbar_ticks=[0,0.25,0.5,0.75,1]
# colmap='Blues'
# cb1=matplotlib.colorbar.ColorbarBase(colorbar, cmap=colmap,ticks=(colbar_ticks)) 
# cb1.set_ticks(colbar_ticks)
# cb1.set_ticklabels(\
# ("{:.2E}".format(Decimal(Min_val)),"{:.3f}".format(0.25*Max_val),\
#  "{:.3f}".format(0.5*Max_val),"{:.3f}".format(0.75*Max_val),\
#  "{:.3f}".format(Max_val)))

# cb1.ax.tick_params(labelsize=size_ticks)
# cb1.outline.set_visible(False)
# cb1.set_label("Relative rate",fontsize=size_axis)
#===============================================

#3.6. (0,1) - Print parameters
#===============================================  
Parameters=plt.subplot(gs[0,1]) 
plt.axis('off') 
plt.text(0, 1, 'r='+"{:.3f}".format(r))
plt.text(0, 0.85, 'K='+"{:.2E}".format(K)) 
plt.text(0, 0.7, 'd='+"{:.2E}".format(d))   
plt.text(0, 0.55, 'c='+"{:.3f}".format(c))
plt.text(0, 0.4, 'm='+"{:.3f}".format(m))  
#===============================================

# #3.7. (2,0) - Plot Tensions
# #===============================================
# ax_20=plt.subplot(gs[2,0])
# plt.xlim([time_0, time_f])
# plt.axhline(y=r/(c*d*B_eq),linewidth=2, color='k',linestyle='--')
# plt.axhline(y=1 -(c*d*B_eq/r) ,linewidth=2, color='k',linestyle='--')
# plt.plot(time, Relative_Tension_Growth,linewidth=3)
# plt.plot(time, Diff_Tension_Growth,linewidth=3)
                                                                               
# plt.legend(['Relative Tension_Eq', r'$\Delta T_{eq}$', 'Relative Tension', r'$\Delta T$'],loc='upper right')
# plt.ylabel('Tensions',fontsize=size_axis)

# #3.4.1. Ticks, labels, and axes
# #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# plt.xticks(fontsize=size_ticks)
# plt.yticks(fontsize=size_ticks)
# #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# #===============================================

# #3.8. (3,0) - Plot Tensions
# #===============================================
# ax_30=plt.subplot(gs[3,0])
# plt.plot(time,Bacteria_Relative_Equilibrium ,linewidth=2,color='orange',label='Bac\teria')
# plt.plot(time,Phage_Relative_Equilibrium ,linewidth=2,color='cyan',label='Phage')
# plt.legend(['Bacteria', 'Phage'], loc='upper right')
# plt.ylabel('% Equilibrium',fontsize=size_axis)
# plt.xlabel(r'time $(h)$',fontsize=size_axis)
# #===============================================

#3.9. Show and save
#===============================================

plt.savefig(path_plot + 'L-V_Hyperion_Test_Sample_'+str(Sampling_Point) + '.pdf',dpi=300)
plt.savefig(path_plot + 'L-V_Hyperion_Test_Sample_'+str(Sampling_Point) + '.png',dpi=300)

plt.show()
#===============================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

