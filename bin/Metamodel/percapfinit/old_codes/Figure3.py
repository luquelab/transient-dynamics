#23/3/2022. Calculate errors for different values of epsilon

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
def Lotka_Volterra_Growth(t,y,Thresholds):
    return [r*y[0], 0]

#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Critical_Concentration_Burst(t,y,Thresholds):
    return y[0] - Thresholds[0]

Critical_Concentration_Burst.terminal=True
Critical_Concentration_Burst.direction = 0
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================

#2. Solver for Lotka-Volterra with growth and burst
#=============================================================
def Lotka_Volterra_Growth_Burst(t,y,Thresholds):    
    return [r*y[0],c*d*y[0]*y[1]]


#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Critical_Concentration_Infection(t,y,Thresholds):
    return y[1] - Thresholds[0]
Critical_Concentration_Infection.terminal=True
Critical_Concentration_Infection.direction = 0
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================

#4. Solver for Lotka-Volterra with growth, infection, and burst
#=============================================================
def Lotka_Volterra_Growth_Infection_Burst(t,y,Thresholds):    
    return [r*y[0] - d*y[0]*y[1], c*d*y[0]*y[1]]

#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Critical_Concentration_Burst_Off(t,y,Thresholds):
    return y[0] - Thresholds[1]
Critical_Concentration_Burst_Off.terminal=True
Critical_Concentration_Burst_Off.direction = 0
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================

#3. Solver for Lotka-Volterra with growth and infection active
#=============================================================
def Lotka_Volterra_Growth_Infection(t,y,Thresholds):    

    return [r*y[0] - d*y[0]*y[1], 0]

#Breaking condition for sensitive bottom
#:::::::::::::::::::::::::::::::::::::::::::::
def Min_Volume_B_Growth_Infection(t,y,Thresholds):

    return y[0]
Min_Volume_B.terminal=True
Min_Volume_B.direction = 1
#:::::::::::::::::::::::::::::::::::::::::::::
#=============================================================

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
def Time_Vector(t_0, t_f,Overlaps=False):

    Resolution=100;Steps=time_f*Resolution + 1
    Step_Size=1/Resolution
    Time_Vector=np.arange(t_0,t_f+Step_Size,Step_Size)
    if Overlaps==True:
        Time_Vector=np.arange(t_0+Step_Size,t_f+Step_Size,Step_Size)

    return Time_Vector
#=============================================================

#5. Time vector for simplified dynamics
#=============================================================
def Solve_Simplified_Model(y0_fun,t_in,Critical_fun,Model_fun,Event_fun):
    step=0.01

    #Solve
    sol_Model_Raw=solve_ivp(lambda t,x:Model_fun(t,x,Critical_fun),[t_in[0],t_in[-1]],y0_fun,method='RK45',dense_output=True,events=lambda t,x: Event_fun(t,x,Critical_fun),max_step=step)

    try:

        #Get events for next stage
        t_f=sol_Model_Raw.t_events[0][0]
        y_f=sol_Model_Raw.y_events[0][0]

        #New time vector
        t_out=Time_Vector(t_in[0],t_f)
        sol_Model=sol_Model_Raw.sol(t_out)
        
    except IndexError:

        t_out=t_in
        sol_Model=sol_Model_Raw.sol(t_out)
        y_f=sol_Model[-1]
    
    return sol_Model, t_out, y_f
#=============================================================


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. MODEL AND SOLVER PARAMETERS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Parameters
#==============================
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

#Full Model
#===================================================================
#2.1. Initial conditions and Events
y0=[1000,10000]

#2.2. Solver
sol_LV=solve_ivp(Lotka_Volterra,[time[0],time[-1]],y0,\
method='RK45',dense_output=True,events=Min_Volume_B,max_step=step)

#2.3. Change time vector in case event triggered
try:
    final_t=sol_LV.t_events[0][0]
    time=time[:int(final_t)]
except IndexError:
    print('No event found')
#===================================================================

#2.4. Calculate individual terms
#====================================================================

#2.4.1. Choose between discrete solution and continuous
z = sol_LV.sol(time) #Discrete solution. Overlaps with time vector
#2.4.2. Initial concentrations
Bacteria_0=y0[0];Phage_0=y0[1]
#2.4.3. Compute rates and extract individual terms
Bioprocesses=PerCapita_Rates_LV(z, time, Parameters)
#2.4.4. Build dataframe with percapita rates, rates, and time
Infection=Bioprocesses['Infection']; Burst=Bioprocesses['Burst']
df=pd.DataFrame({'Time':time,'Bacterial_Concentration':z[0],'Phage_Concentration':z[1],'Per_Capita_Infection':Infection,'Per_Capita_Burst':Burst})

#2.4.5. Find critical concentrations and times
#.............................................................................
#Define epsilon
epsilon_1=0.1;epsilon_2=1
epsilon=epsilon_2

epsilon_vector=np.arange(0.1,1.4,0.1)
print(type(epsilon_vector))
epsilon_vector=[0.005,0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5,0.75, 1]
# epsilon_vector=np.asarray(epsilon_vector)
# print(type(epsilon_vector))
#print(np.logspace(-3, 0, 6, endpoint=False))

Mean_Relative_Errors_Bacteria=[]
Mean_Relative_Errors_Phage=[]
Mean_Relative_Errors=[]

for epsilon_i in epsilon_vector:
    print(epsilon_i)
    #Find critical indices for burst
    Per_Capita_Burst=df['Per_Capita_Burst'].to_numpy()
    Critical_Indices_Burst=Epsilon_Finder(Per_Capita_Burst, epsilon_i)
    #Find critical concentrations and times
    Critical_Concentrations_Burst=[];Critical_Times_Burst=[]
    for index in Critical_Indices_Burst:
        Critical_Concentrations_Burst.append(df.iloc[index]['Bacterial_Concentration'])

    #Find critical indices for infection
    Per_Capita_Infection=df['Per_Capita_Infection'].to_numpy()
    Critical_Indices_Infection=Epsilon_Finder(Per_Capita_Infection, epsilon_i)
    print(Critical_Indices_Infection)

    #Find critical concentration
    Critical_Concentrations_Infection=df.iloc[int(Critical_Indices_Infection)]['Phage_Concentration']
    Critical_Concentrations_Infection=[Critical_Concentrations_Infection]
    #Find critical times
    # Critical_Infection_Time=df.iloc[int(Critical_Indices_Infection)]['Time']
    # print(Critical_Infection_Time)
    #.............................................................................
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Solve system with activated/deactivated terms
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #1. Only growth 
    #.............................................................
    #1.1. Initial conditions
    y0=[1000,10000]

    #Solution, time, and events
    z_g, time_g, y_f_growth=Solve_Simplified_Model(y0, time, Critical_Concentrations_Burst, Lotka_Volterra_Growth, Critical_Concentration_Burst)
    #.............................................................

    #2. Growth + Burst
    #.............................................................
    #2.1. Time vector and initial conditions
    time_0_g_b=time_g[-1];time_g_b=Time_Vector(time_0_g_b,time_f,Overlaps=True)
    y0_g_b=[y_f_growth[0],y_f_growth[1]]

    #Solution, time, and events
    z_g_b, time_g_b, y_f_g_b=Solve_Simplified_Model(y0_g_b, time_g_b, Critical_Concentrations_Infection, Lotka_Volterra_Growth_Burst, Critical_Concentration_Infection)
    #.............................................................

    #3. Growth, Infection, and burst
    #.............................................................
    #3.1.Time vector and initial conditions
    time_0_g_b_i=time_g_b[-1];
    time_g_b_i=Time_Vector(time_0_g_b_i,time_f,Overlaps=True)
    
    y0_g_b_i=[y_f_g_b[0],y_f_g_b[1]]

    #Solution, time, and events
    z_g_b_i,time_g_b_i,y_f_g_b_i=Solve_Simplified_Model(y0_g_b_i, time_g_b_i, Critical_Concentrations_Burst, Lotka_Volterra_Growth_Infection_Burst, Critical_Concentration_Burst_Off)
    #.............................................................

    #4. Growth, Infection
    #.............................................................
    #4.1. Time vector and initial conditions
    time_0_g_i=time_g_b_i[-1];
    time_g_i=Time_Vector(time_0_g_i,time_f,Overlaps=True)
    y0_g_i=[y_f_g_b_i[0],y_f_g_b_i[1]]

    #Solution, time, and events
    z_g_i,time_g_i,y_f_g_b_i=Solve_Simplified_Model(y0_g_i, time_g_i, Critical_Concentrations_Burst, Lotka_Volterra_Growth_Infection,Min_Volume_B_Growth_Infection)
    #.............................................................

    Simplified_Dynamics=np.concatenate((z_g,z_g_b,z_g_b_i,z_g_i), axis=1)
    Time_Simplified_Dynamics=np.concatenate((time_g,time_g_b,time_g_b_i,time_g_i))

    #Relative error bacteria
    #---------------------------------------------------------------
    Relative_Error_Bacteria=[]
    for (bacteria_true, bacteria_model) in zip(z[0], Simplified_Dynamics[0]):
        Relative_Error_Bacteria_i=\
        (np.absolute(bacteria_true-bacteria_model))/bacteria_true 
        Relative_Error_Bacteria.append(Relative_Error_Bacteria_i)

    Mean_Relative_Error_Bacteria=np.mean(Relative_Error_Bacteria)
    #---------------------------------------------------------------

    #Relative error phage
    #---------------------------------------------------------------
    Relative_Error_Phage=[]
    for (bacteria_true, bacteria_model) in zip(z[1], Simplified_Dynamics[1]): 
        Relative_Error_Phage_i=\
        (np.absolute(bacteria_true-bacteria_model))/bacteria_true
        Relative_Error_Phage.append(Relative_Error_Phage_i)

    Mean_Relative_Error_Phage=np.mean(Relative_Error_Phage)
    #---------------------------------------------------------------
    
    #Total mean Error
    #---------------------------------------------------------------
    Mean_Relative_Error=np.mean([Mean_Relative_Error_Phage,Mean_Relative_Error_Bacteria])

    print("Error")
    print(Mean_Relative_Error)
    #---------------------------------------------------------------

    Mean_Relative_Errors_Bacteria.append(Mean_Relative_Error_Bacteria)
    Mean_Relative_Errors_Phage.append(Mean_Relative_Error_Phage)
    Mean_Relative_Errors.append(Mean_Relative_Error)

    
#Compute relative error by sections 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Burst_on=len(time_g)
    Infect_on=len(time_g) + len(time_g_b)
    Burst_off=len(time_g) + len(time_g_b) +len(time_g_b_i)
    End=len(time_g) + len(time_g_b) +len(time_g_b_i) + len(time_g_i)
    
    #Relative error phage
    #--------------------------------------------------------------------------------------------
    Tramos_Dynamics_Model_P=[z_g[1],z_g_b[1],z_g_b_i[1],z_g_i[1]]
    Tramos_Dynamics_True_P=[z[1][:Burst_on],z[1][Burst_on:Infect_on],z[1][Infect_on:Burst_off],z[1][Burst_off:]]

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
    #--------------------------------------------------------------------------------------------

    #Relative error bacteria
    #--------------------------------------------------------------------------------------------
    Tramos_Dynamics_Model_B=[z_g[0],z_g_b[0],z_g_b_i[0],z_g_i[0]]
    Tramos_Dynamics_True_B=[z[0][:Burst_on],z[0][Burst_on:Infect_on],z[0][Infect_on:Burst_off],z[0][Burst_off:]]
    
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
    #--------------------------------------------------------------------------------------------
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Plot Figures
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Path to save figure
Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/Paper_Matt/'

plt.plot(time,Relative_Error_Bact_Total,color='b',linewidth=1,label='Total Error')
plt.plot(time,Relative_Error_Phage_Total,color='r',linewidth=1,label='Total Error')
plt.axvline(time_g[-1],color='k',linewidth=1,linestyle='dotted')
plt.axvline(time_g_b[-1],color='k',linewidth=1,linestyle='dotted')
plt.axvline(time_g_b_i[-1],color='k',linewidth=1,linestyle='dotted')
plt.show()

#Fontsizes
size_axis=7;size_ticks=6;size_title=5

#Figure Size
cm = 1/2.54  # centimeters in inches
Width=8*cm;Height=4*cm #Width and height of plots
figure(figsize=(Width, Height), dpi=300)

#Axes, title and ticks
plt.ylabel('Mean Relative Error',fontsize=size_axis)
plt.xlabel(r'$\epsilon$ threshold',fontsize=size_axis)

plt.plot(epsilon_vector,Mean_Relative_Errors,color='k',linewidth=1,label='Total Error')
plt.plot(epsilon_vector,Mean_Relative_Errors_Bacteria,color='b',linewidth=1, label='Bacterial Error')
plt.plot(epsilon_vector,Mean_Relative_Errors_Phage,color='r',linewidth=1, label='Phage Error')
plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)
plt.legend(loc='upper left',fontsize=size_ticks, frameon=False)

plt.xscale('log')
sns.despine(top=True, right=True, left=False, bottom=False)

plt.show()


