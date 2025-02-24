#24/06/2021 - This function solves the Lotka-Volterra system

# Import libraries
#++++++++++++++++++++++++++++++++
import numpy as np
from scipy.integrate import odeint
import math
from scipy.integrate import solve_ivp
import pandas as pd
from decimal import Decimal
#++++++++++++++++++++++++++++++++

#0. FUNCTIONS.
#0.0. Solution to Lotka-Volterra
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Solve_Lotka_Volterra(Parameters, Initial_Conditions, time_vector, step):
    
    r=Parameters['r']
    K=Parameters['K']
    d=Parameters['d']
    c=Parameters['c']
    m=Parameters['m']
    
    #0.1. Lotka-Volterra system
    #=================================================    
    def Lotka_Volterra(t,y):

        return [r*(1-(y[0]/float(K)))*y[0] - d*y[0]*y[1],
                c*d*y[0]*y[1] - m*y[1]]

    #Breaking conditions for sensitive bacteria. Bottom
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #Critical Concentrations
    #----------------------------------------
    def Critical_Concentration_B_growth(t,y):
        return y[0]-K
    Critical_Concentration_B_growth.direction = 0

    def Critical_Concentration_P_infection(t,y):
        return y[1]-r/d
    Critical_Concentration_P_infection.direction = 0

    def Critical_Concentration_B_burst(t,y):
        return y[0]-m/c*d
    Critical_Concentration_B_burst.direction = 0
    #----------------------------------------
    
    def Min_Volume_B(t,y):
        return y[0] - 1
    Min_Volume_B.terminal=True
    Min_Volume_B.direction = 1
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #=================================================

    #0.2. MODEL SOLUTION
    #====================================================================
    Events_Model=[Critical_Concentration_B_growth, \
    Critical_Concentration_P_infection,Critical_Concentration_B_burst]
    
    y0=[Initial_Conditions[0], Initial_Conditions[1] ]

    #0.2.1. Solver
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    sol_LV=solve_ivp(Lotka_Volterra,[time_vector[0],time_vector[-1]],y0,\
    method='RK45',dense_output=True,events=Events_Model,max_step=step)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.2.2. Change time vector in case event triggered
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    try:
        final_t=sol_LV.t_events[0][0]
        time=time_vector[:int(time_vector[-1])]
        print('Event found')
    except IndexError:
        pass
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    return sol_LV
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#0.1. Solution to Experimental System
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Solve_Experiment_Induction(Parameters, Initial_Conditions, time_vector, step):
    
    r=Parameters['r']
    K=Parameters['K']
#    d=Parameters['d']
    c=Parameters['c']
    m=Parameters['m']
    mu=Parameters['mu']
    
    #0.1. Induction system
    #=================================================    
    def Lotka_Volterra(t,y):

        return [r*(1-(y[0]/float(K)))*y[0] - mu*y[0],
                c*mu*y[0] - m*y[1]]

    #Breaking conditions for sensitive bacteria. Bottom
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #Critical Concentrations
    #----------------------------------------
    def Critical_Concentration_L_growth(t,y):
        return y[0]-K
    Critical_Concentration_L_growth.direction = 0

    # def Critical_Concentration_P_infection(t,y):
    #     return y[1]-r/d
    # Critical_Concentration_P_infection.direction = 0

    # def Critical_Concentration_B_burst(t,y):
    #     return y[0]-m/c*d
    # Critical_Concentration_B_burst.direction = 0
    #----------------------------------------
    
    def Min_Volume_L(t,y):
        return y[0] - 1
    Min_Volume_B.terminal=True
    Min_Volume_B.direction = 1
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #=================================================

    #0.2. MODEL SOLUTION
    #====================================================================
    Events_Model=[Critical_Concentration_B_growth, \
    Critical_Concentration_P_infection,Critical_Concentration_B_burst]
    
    y0=[Initial_Conditions[0], Initial_Conditions[1] ]

    #0.2.1. Solver
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    sol_Exp=solve_ivp(Experiment_Induction,[time_vector[0],time_vector[-1]],y0,\
    method='RK45',dense_output=True,events=Events_Model,max_step=step)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.2.2. Change time vector in case event triggered
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    try:
        final_t=sol_Exp.t_events[0][0]
        time=time_vector[:int(time_vector[-1])]
        print('Event found')
    except IndexError:
        pass
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    return sol_Exp
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#0.1. Split range of parameters in equispaced or logspaced intervals
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def LHS_intervals(min_value,max_value,Samples_fun,Logarithmic=0):
    interval=np.linspace(min_value, max_value, Samples_fun+1)
    
    if Logarithmic==1:
        interval=np.logspace(np.log10(min_value), np.log10(max_value),\
                             num=Samples_fun+1, endpoint=True, base=10)
        
    return interval 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.2. Latin_Hypercube_Sampling
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
def Latin_Hypercube_Sampling(Ranges_Parameters, Sampling_Points, Seed):
    np.random.seed(Seed)

    #A dictionary that gives an id from 1 to N to each interval
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Ids_Intervals={}  
    for Parameter in Ranges_Parameters:
        
        Min_Value=Ranges_Parameters[Parameter][0]
        Max_Value=Ranges_Parameters[Parameter][1]
        vector_interval=LHS_intervals(Min_Value, Max_Value, Sampling_Points)

        #Condition for logarithmic sampling
        #---------------------------------------------------------------------
        Min_o_magnitude=np.log10(Min_Value)     
        Max_o_magnitude=np.log10(Max_Value)
        
        Delta_Orders_of_Magnitude=\
                    abs(abs(np.log10(Max_Value)) - abs(np.log10(Min_Value)))
        
        if Delta_Orders_of_Magnitude>=2:
            vector_interval=LHS_intervals(Min_Value, Max_Value, Sampling_Points,Logarithmic=1) 
        #---------------------------------------------------------------------

        for i in range(len(vector_interval)-1):
            try:
                Ids_Intervals[Parameter][i]=\
                                    [vector_interval[i],vector_interval[i+1]]
            except KeyError:
                Ids_Intervals[Parameter]=\
                                    {i:[vector_interval[i],vector_interval[i+1]]}
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Randomly choose intervals for the parameters
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Samples={}
    for parameter in Ranges_Parameters:
        Samples[parameter]=\
        list(np.random.choice(Sampling_Points,Sampling_Points,replace=False))
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Uniformly choose random values within intervals and store them in dictionary
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    counter_samples=0
    Latin_Hypercube_Sampling={}

    for parameter in Samples:
        for interval in Samples[parameter]:
            #Interval minimum
            Min_value=Ids_Intervals[parameter][interval][0]
            #Interval maximum
            Max_value=Ids_Intervals[parameter][interval][1]
            #Uniformly choose random value
            parameter_value=np.random.uniform(Min_value,Max_value)

            #Store value in dictionary
            #---------------------------------------------------------------------
            try:
                Latin_Hypercube_Sampling[counter_samples][parameter]=\
                                                    parameter_value
            except KeyError:
                Latin_Hypercube_Sampling[counter_samples]=\
                                            {parameter:parameter_value}
            #--------------------------------------------------------------------
            counter_samples+=1
            if counter_samples==Sampling_Points:
                counter_samples=0
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    return Latin_Hypercube_Sampling
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
