#14/03/2022. Multiple values of epsilon to compute the error between the actual solution and the simplified solution

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
#++++++++++++++++++++++++++++++++++++++++++++

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0. Solver for Lotka-Volterra
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

#1. Solver for Lotka-Volterra with only growth active
#=============================================================
def Lotka_Volterra_Growth(t,y):    

    return [r*y[0], 0]


#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Min_Volume_B_Growth(t,y):

    return y[0] - 1
Min_Volume_B.terminal=True
Min_Volume_B.direction = 1
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================

#2. Solver for Lotka-Volterra with growth and burst
#=============================================================
def Lotka_Volterra_Growth_Burst(t,y):    

    return [r*y[0],c*d*y[0]*y[1]]


#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Min_Volume_B_Growth_Infection(t,y):

    return y[0] - 1
Min_Volume_B.terminal=True
Min_Volume_B.direction = 1
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================


#2. Solver for Lotka-Volterra with only growth active
#=============================================================
def Lotka_Volterra_Growth_Infection(t,y):    

    return [r*y[0] - d*y[0]*y[1], 0]


#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Min_Volume_B_Growth_Infection(t,y):

    return y[0] - 1
Min_Volume_B.terminal=True
Min_Volume_B.direction = 1
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================


#3. Solver for Lotka-Volterra with growth, infection, and burst
#=============================================================
def Lotka_Volterra_Growth_Infection_Burst(t,y):    

    return [r*y[0] - d*y[0]*y[1], c*d*y[0]*y[1]]


#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Min_Volume_B_Growth_Infection_Burst(t,y):

    return y[0] - 1
Min_Volume_B.terminal=True
Min_Volume_B.direction = 1
#:::::::::::::::::::::::::::::::::::::::::::::


#4. Find thresholds for epsilon
#=============================================================
def Epsilon_Finder(Per_Capita_Process, epsilon_value):
    Critical_Indices_Per_Capita=\
    np.where(np.sign(Per_Capita_Process[:-1] - epsilon_value) != \
    np.sign(Per_Capita_Process[1:] - epsilon_value) )[0] + 1

    return Critical_Indices_Per_Capita
#=============================================================


#5. Time vector for simplified dynamics
#=============================================================
def Time_Vector(t_0, t_f):

    Resolution=100;Steps=time_f*Resolution + 1
    Step_Size=1/Resolution
    Time_Vector=np.arange(t_0,t_f+Step_Size,Step_Size)

    return Time_Vector
#=============================================================

#1. MODEL AND SOLVER PARAMETERS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Initial conditions
#==============================
S0_values=[1e3, 1e6, 1e9]
T0_values=[1e3, 1e6, 1e9]
#==============================

#1.2. Parameters
#==============================

#Marine Parameters
#::::::::::::::::::::
# r=1/float(24)
# K=1e7
# d=1e-8
# c=5
# m=1/float(6)
#m=r
#::::::::::::::::::::

#T4-E. coli parameters
#::::::::::::::::::::
r=0.9
K=1e7
d=3e-8
c=150
m=2.8e-3
Volume=1
Parameters={'r':r,'d':d,'c':c,'m':m,'K':K}
#::::::::::::::::::::
#==============================

#1.3. Solver Parameters
#==============================
time_0=0
time_f=14
Resolution=100
Steps=time_f*Resolution + 1
time=np.linspace(time_0,time_f,Steps)
step=0.01
#==============================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#2. MODEL SOLUTION
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.1. Initial conditions and Events
#==============================
Events_Model=[Min_Volume_B]
y0=[1000,10000]
#==============================

#2.2. Solver
#====================================================================
sol_LV=solve_ivp(Lotka_Volterra,[time[0],time[-1]],y0,\
method='RK45',dense_output=True,events=Events_Model,max_step=step)
#====================================================================

#2.3. Change time vector in case event triggered
#=====================================
try:
    final_t=sol_LV.t_events[0][0]
    time=time[:int(final_t)]
except IndexError:
    print('No event found')
#=====================================


#2.4. Calculate individual terms
#====================================================================

#2.4.1. Choose between discrete solution and continuous
#:::::::::::::::::::::::::::::::::
z = sol_LV.sol(time) #Discrete solution. Overlaps with time vector
#:::::::::::::::::::::::::::::::::

#2.4.2. Initial concentrations
#::::::::::::::::::::::::::::::
Bacteria_0=y0[0]
Phage_0=y0[1]
#::::::::::::::::::::::::::::::

#2.4.3. Compute rates and extract individual terms
Bioprocesses=PerCapita_Rates_LV(z, time, Parameters)

#2.4.4. Build dataframe with percapita rates, rates, and time
Infection=Bioprocesses['Infection']; Burst=Bioprocesses['Burst']

df=pd.DataFrame({'Time':time,'Bacterial_Concentration':z[0],'Phage_Concentration':z[1],'Per_Capita_Infection':Infection,'Per_Capita_Burst':Burst})

#2.4.5. Find critical concentrations and times
#.............................................................................
#Define epsilon
epsilon_1=0.1;epsilon_2=1
epsilon=epsilon_1
step=0.01

epsilon_vector=np.arange(0.1,1.4,0.1)
Mean_Relative_Errors=[];Mean_Relative_Errors_Bacteria=[];
Mean_Relative_Errors_Phage=[]
for epsilon_i in epsilon_vector:
    #Find critical indices for burst
    
    Per_Capita_Burst=df['Per_Capita_Burst'].to_numpy()
    Critical_Indices_Burst=Epsilon_Finder(Per_Capita_Burst, epsilon_i)
#    print(Critical_Indices_Burst)

    Critical_Concentrations_Burst=[]; Critical_Times_Burst=[]
    for index in Critical_Indices_Burst:
        Critical_Concentrations_Burst.append(df.iloc[index]['Bacterial_Concentration'])
        Critical_Times_Burst.append(df.iloc[index]['Time'])

    Activation_Burst_Time=Critical_Times_Burst[0]
    Inactivation_Burst_Time=Critical_Times_Burst[1]

    #Find critical indices for infection
    Per_Capita_Infection=df['Per_Capita_Infection'].to_numpy()
    Critical_Indices_Infection=Epsilon_Finder(Per_Capita_Infection, epsilon_i)
    #Find critical times
    Critical_Infection_Time=df.iloc[int(Critical_Indices_Infection)]['Time']
    print(Critical_Infection_Time)    
    #.............................................................................
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Solve system with activated/deactivated terms
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    #1. Only growth 
    #.............................................................
    #Time_g
    #==============================
    time_0_g=0.0;time_f_g=Activation_Burst_Time
    time_g=Time_Vector(time_0_g, Activation_Burst_Time)
    #==============================

    #2. Initial conditions and Events
    #==============================
    Events_Model=[Min_Volume_B]
    y0=[1000,10000]
    #==============================

    #3. Solver
    #====================================================================
    sol_LV_G=solve_ivp(Lotka_Volterra_Growth,[time_g[0],time_g[-1]],y0,\
    method='RK45',dense_output=True,events=Events_Model,max_step=step)
    z_g=sol_LV_G.sol(time_g) #Discrete solution. Overlaps with time vector
    #====================================================================

    #2. Growth + Burst
    #.............................................................
    #Time_g_b
    #==============================
    time_0_g_b=Activation_Burst_Time;time_f_g_b=Critical_Infection_Time
    time_g_b=Time_Vector(time_0_g_b,time_f_g_b)
    #==============================

    #2. Initial conditions and Events
    #===============================
    Events_Model=[Min_Volume_B]
    Bacteria_0_g_b=z_g[0][-1]
    Phage_0_g_b=z_g[1][-1]
    y0_g_b=[Bacteria_0_g_b,Phage_0_g_b]
    #===============================

    #3. Solver
    #====================================================================
    sol_LV_G_B=solve_ivp(Lotka_Volterra_Growth_Burst,[time_g_b[0],time_g_b[-1]],y0_g_b,method='RK45',dense_output=True,events=Events_Model,max_step=step)

    z_g_b=sol_LV_G_B.sol(time_g_b) #Discrete solution. Overlaps with time vector
    #====================================================================
    #.............................................................

    #3. Growth, Infection, and burst
    #.............................................................
    #Time
    #==============================
    time_0_g_b_i=Critical_Infection_Time;time_f_g_b_i=Inactivation_Burst_Time
    time_g_b_i=Time_Vector(time_0_g_b_i,time_f_g_b_i)
    #==============================
    
    #2. Initial conditions and Events
    #===============================
    Events_Model=[Min_Volume_B]
    Bacteria_0_g_b_i=z_g_b[0][-1]
    Phage_0_g_b_i=z_g_b[1][-1]
    y0_g_b_i=[Bacteria_0_g_b_i,Phage_0_g_b_i]
    #===============================

    #3. Solver
    #====================================================================
    sol_LV_G_B_I=solve_ivp(Lotka_Volterra_Growth_Infection_Burst,[time_g_b_i[0],time_g_b_i[-1]],y0_g_b_i,method='RK45',dense_output=True,events=Events_Model,max_step=step)

    z_g_b_i=sol_LV_G_B_I.sol(time_g_b_i) #Discrete solution. Overlaps with time vector
    #====================================================================
    #.............................................................

    #4. Growth, Infection
    #.............................................................
    #Time
    #==============================
    time_0_g_i=Inactivation_Burst_Time;time_f_g_i=14
    time_g_i=Time_Vector(time_0_g_i,time_f_g_i)
    #==============================

    #2. Initial conditions and Events
    #===============================
    Events_Model=[Min_Volume_B]
    Bacteria_0_g_i=z_g_b_i[0][-1]
    Phage_0_g_i=z_g_b_i[1][-1]
    y0_g_i=[Bacteria_0_g_i,Phage_0_g_i]
    #===============================

    #3. Solver
    #====================================================================
    sol_LV_G_I=solve_ivp(Lotka_Volterra_Growth_Infection,[time_g_i[0],time_g_i[-1]],y0_g_i,method='RK45',dense_output=True,events=Events_Model,max_step=step)
    z_g_i=sol_LV_G_I.sol(time_g_i) #Discrete solution. Overlaps with time vector
    #====================================================================


    Simplified_Dynamics=np.concatenate((z_g,z_g_b,z_g_b_i,z_g_i), axis=1)
    Time_Simplified_Dynamics=np.concatenate((time_g,time_g_b,time_g_b_i,time_g_i))

    Relative_Error_Bacteria=[]
    for (bacteria_true, bacteria_model) in zip(z[0], Simplified_Dynamics[0]):
        Relative_Error_Bacteria_i=(np.absolute(bacteria_true - bacteria_model))/bacteria_true
        Relative_Error_Bacteria.append(Relative_Error_Bacteria_i)
    Mean_Relative_Error_Bacteria=np.mean(Relative_Error_Bacteria)


    Relative_Error_Phage=[]
    for (bacteria_true, bacteria_model) in zip(z[1], Simplified_Dynamics[1]):
        Relative_Error_Phage_i=(np.absolute(bacteria_true - bacteria_model))/bacteria_true
        Relative_Error_Phage.append(Relative_Error_Phage_i)
    Mean_Relative_Error_Phage=np.mean(Relative_Error_Phage)

    Mean_Relative_Error=np.mean([Mean_Relative_Error_Phage,Mean_Relative_Error_Bacteria])
#    print(Mean_Relative_Error)
    Mean_Relative_Errors.append(Mean_Relative_Error)
    Mean_Relative_Errors_Phage.append(Mean_Relative_Error_Phage)
    Mean_Relative_Errors_Bacteria.append(Mean_Relative_Error_Bacteria)



    #Compute relative error by sections
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    corte1=len(time_g)
    corte2=len(time_g) + len(time_g_b)
    corte3=len(time_g) + len(time_g_b) +len(time_g_b_i)
    corte4=len(time_g) + len(time_g_b) +len(time_g_b_i) + len(time_g_i)

    Tramos_Dynamics_Model_P=[z_g[1],z_g_b[1],z_g_b_i[1],z_g_i[1]] 
    Tramos_Dynamics_True_P=[z[1][:corte1],z[1][corte1:corte2],z[1][corte2:corte3],z[1][corte3:]]

    Relative_Error_Phage_Total=[]#Error of the whole vector
    Relative_Error_Phage_Vec=[]#Vector of errors of the whole vector
    for (din1, din2) in zip(Tramos_Dynamics_True_P,Tramos_Dynamics_Model_P):
        Relative_Error_Phage_Tramo=[]
        for (bacteria_true, bacteria_model) in zip(din1, din2): 
            Relative_Error_Phage_i=\
            (np.absolute(bacteria_true-bacteria_model))/bacteria_true             
            Relative_Error_Phage_Tramo.append(Relative_Error_Phage_i)
            Relative_Error_Phage_Total.append(Relative_Error_Phage_i)
        Relative_Error_Phage_Vec.append(sum(Relative_Error_Phage_Tramo)/len(z[0]))


    Tramos_Dynamics_Model_B=[z_g[0],z_g_b[0],z_g_b_i[0],z_g_i[0]] 
    Tramos_Dynamics_True_B=[z[0][:corte1],z[0][corte1:corte2],z[0][corte2:corte3],z[0][corte3:]]
    
    Relative_Error_Bact_Total=[]#Error of the whole vector
    Relative_Error_Bact_Vec=[]#Vector of errors of the whole vector
    for (din1, din2) in zip(Tramos_Dynamics_True_B,Tramos_Dynamics_Model_B):
        Relative_Error_Bact_Tramo=[]
        for (bacteria_true, bacteria_model) in zip(din1, din2):
            Relative_Error_Bact_i=\
            (np.absolute(bacteria_true-bacteria_model))/bacteria_true
            Relative_Error_Bact_Tramo.append(Relative_Error_Bact_i)
            Relative_Error_Bact_Total.append(Relative_Error_Bact_i)

        Relative_Error_Bact_Vec.append(sum(Relative_Error_Bact_Tramo)/len(z[0]))
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

print(Relative_Error_Bact_Vec)
print(Relative_Error_Phage_Vec)



#Fontsizes 
size_axis=7;size_ticks=5;size_title=5 
plt.plot(epsilon_vector, Mean_Relative_Errors,color='k',linewidth=2)
#plt.plot(epsilon_vector, Mean_Relative_Errors_Bacteria,color='b',linewidth=2)
#plt.plot(epsilon_vector, Mean_Relative_Errors_Phage,color='r',linewidth=2)

#Axes, title and ticks
plt.ylabel('Mean Relative Error',fontsize=size_axis)
plt.xlabel(r'$\epsilon$',fontsize=size_axis)
plt.show()



#Plot Figures
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Path to save figure
Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/Paper_Matt/'

plt.plot(time,Relative_Error_Bact_Total,color='b',linewidth=1,label='Total Error')
plt.plot(time,Relative_Error_Phage_Total,color='r',linewidth=1,label='Total Error'
)                                                                                 
plt.axvline(time_g[-1],color='k',linewidth=1,linestyle='dotted')
plt.axvline(time_g_b[-1],color='k',linewidth=1,linestyle='dotted')
plt.axvline(time_g_b_i[-1],color='k',linewidth=1,linestyle='dotted')
plt.show()
#Fontsizes
size_axis=7;size_ticks=5;size_title=5
#Figure Size
cm = 1/2.54  # centimeters in inches
Width=8*cm;Height=4*cm #Width and height of plots
figure(figsize=(Width, Height), dpi=300)
#Colors
cmap='RdBu';cmap_pieces= matplotlib.cm.get_cmap(cmap)
Phage_color=cmap_pieces(0.1);Bacteria_color=cmap_pieces(0.9)

#Plot Figure
plt.plot(time,z[0].T,color=Bacteria_color,linewidth=1,label='Prey - Full model')
plt.plot(time,z[1].T,color=Phage_color,linewidth=1,label='Predator - Full model')

#Growth
plt.plot(time_g, z_g[0].T,color='k',linestyle='--',dashes=(5, 5),linewidth=1)
plt.plot(time_g, z_g[1].T,color='k',linestyle='--',dashes=(5, 5),linewidth=1)
#Growth+Burst
plt.plot(time_g_b, z_g_b[0].T,color='k',linestyle='--',dashes=(5, 5),linewidth=1)
plt.plot(time_g_b, z_g_b[1].T,color='k',linestyle='--',dashes=(5, 5),linewidth=1)
#Growth+Burst+Infection
plt.plot(time_g_b_i, z_g_b_i[0].T,color='k',linestyle='--',dashes=(5, 5),linewidth=1)
plt.plot(time_g_b_i, z_g_b_i[1].T,color='k',linestyle='--',dashes=(5, 5),linewidth=1)
#Growth+Infection
plt.plot(time_g_i, z_g_i[0].T,color='k',linestyle='--',dashes=(5, 5),linewidth=1,label='Prey - Simplified model')
plt.plot(time_g_i, z_g_i[1].T,color='k',linestyle='--',dashes=(5, 5),linewidth=1,label='Predator - Simplified model')

#Vertical lines indicating critical times
plt.axvline(Critical_Infection_Time,color='k',linewidth=1,linestyle='dotted')
plt.axvline(Inactivation_Burst_Time,color='k',linewidth=1,linestyle='dotted')
plt.axvline(Activation_Burst_Time,color='k',linewidth=1,linestyle='dotted')

print(Critical_Infection_Time)
print(Inactivation_Burst_Time)
print(Activation_Burst_Time)

#Legend and scale
plt.yscale('log')
plt.ylim([1,2e10])
plt.xlim([0,14])

#Axes, title and ticks
plt.ylabel('Concentration $(ml^{-1})$',fontsize=size_axis)
plt.xlabel('time (h)',fontsize=size_axis)
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)
plt.legend(loc='upper left',fontsize=size_ticks)

#Ticks
Time_Dataframe=df.index.values.tolist()
xticks_pos=[position for position in range(Time_Dataframe[0],Time_Dataframe[-1]+1,2*Resolution)]
#print(xticks_pos)
xticks_labels=[label for label in range(int(time_0),int(time_f+1),2)]
#print(xticks_labels)
plt.xticks(xticks_labels)
plt.xticks(rotation='horizontal',fontsize=size_ticks)
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks) 

#Save and show plots
Name_of_Figure='Fig2b_epsilon='+str(epsilon)
Extensions=['.pdf','.eps','.png','.svg']
for ext in Extensions:
    plt.savefig(Output_Path+Name_of_Figure + ext,dpi=300)
plt.show()
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
