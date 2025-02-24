#24/06/2021 - This file lives in the folder /home/sergio/.local/lib/python3.6/site-packages/
#This is a small library with functions that have been useful for my project

# Import libraries
#++++++++++++++++++++++++++++++++
import numpy as np
from scipy.integrate import odeint
import math
from scipy.integrate import solve_ivp
import pandas as pd
from decimal import Decimal
import matplotlib.pyplot as plt
import seaborn as sns
#++++++++++++++++++++++++++++++++

#Index. Functions in this file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#0. TRANSIENT DYNAMICS FUNCTIONS.

#Simplified versions of Lotka-Volterra
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. Growth Active
#=============================================================                                      
def Lotka_Volterra_G(t,y,Thresholds):                                                               
    return [r*y[0],0]                                                                               
#============================================================= 


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#0.1. Extract per capita rates from Lotka-Volterra
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def PerCapita_Rates_Timescale(Concentrations, Time, Parameters, timescale_var):
    
    r=Parameters['r'];K=Parameters['K'];c=Parameters['c']
    d=Parameters['d'];m=Parameters['m']

    Bioprocesses={'Growth':[],'Infection':[],'Burst':[],'Decay':[]}

    for (Bacteria,Phage) in zip(Concentrations[0],Concentrations[1]):

        #Compute global and per capita rates
        #............................................................
        Growth_t=(r*(Bacteria))
        Normalized_Growth_t=timescale_var*Growth_t/(Bacteria)
        
        Infection_t=(d*Bacteria*Phage)
        Normalized_Infection_t=timescale_var*d*Phage
#        Normalized_Infection_t=timescale_var*Infection_t/(Bacteria)
        
        
        Burst_t=(c*d*Bacteria*Phage)
        Normalized_Burst_t=timescale_var*Burst_t/(Phage)
        
        Decay_t=(m*Phage)
        Normalized_Decay_t=timescale_var*Decay_t/(Phage)
        #............................................................

        #Append per capita rates to dictionary
        #............................................................
        Bioprocesses['Growth'].append(Normalized_Growth_t)
        Bioprocesses['Infection'].append(Normalized_Infection_t)
        Bioprocesses['Burst'].append(Normalized_Burst_t)
        Bioprocesses['Decay'].append(Normalized_Decay_t)
        #............................................................

        #Update previous-step concentrations for rates
        #.......................
        Bacteria_0=Bacteria
        Phage_0=Phage
        #.......................

    return Bioprocesses
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Obtain hash for the initial dynamics
#=============================================================
def Get_hash_of_initial_dynamics(Dynamics,Terms_Ordered_fun):
    Active_Terms_List=[]  
    for term in Terms_Ordered_fun:
        Active_Terms_List.append(Dynamics[term])                     
        
    Initial_Hash=''.join(str(value) for value in Active_Terms_List)  
    
    return Initial_Hash                                              
#============================================================= 

#Monitor dynamics
#=============================================================
def Monitor_Dynamics(Dynamics_Dict,Index_var, Terms_Ordered_var):
    
    if Index_var==0:                                                                                
        Dynamics_Dict['Predation']=Dynamics_Dict['Predation']-1                                     
        Dynamics_Dict['Predation']=abs(Dynamics_Dict['Predation'])                                  
        
    elif Index_var==1:                                                                              
        Dynamics_Dict['Burst']=Dynamics_Dict['Burst']-1                                             
        Dynamics_Dict['Burst']=abs(Dynamics_Dict['Burst'])
        
    else:                                                                                           
        Dynamics_Dict=Dynamics_Dict                                                                 
        
    Hash_Dynamics_0=[]                                                                              
    for term in Terms_Ordered_var: 
        Hash_Dynamics_0.append(Dynamics_Dict[term])                                                 
        Hash_Dynamics=''.join(str(value) for value in Hash_Dynamics_0)                              
        
    return Hash_Dynamics, Dynamics_Dict
#============================================================= 

#Build time vector
#=============================================================
def Build_Time_Vector(t_0, t_f):                                                                    
    
    Resolution=100;Steps=t_f*Resolution + 1 
    Step_Size=1/Resolution
    Time_Vector=np.arange(t_0,t_f,Step_Size)
    
    return Time_Vector
#============================================================= 

#Solve Simplified Model
#=============================================================
def Solve_Simplified_Model(y0_fun0,t0_fun,tf_fun,Event_thresholds,Model_fun,Event_functions):
    
    #Build time vector and set initial conditions                                                   
    step=0.01                                                                                       
    t_in=Build_Time_Vector(t0_fun,tf_fun)                                                           
    y0_fun=[ y0_fun0[0], y0_fun0[1] ]                                                               
    
    #Solve simplified model                                                                         
    sol_Model_Raw=solve_ivp(lambda t,x:Model_fun(t,x,Event_thresholds),[t_in[0],t_in[-1]],\
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

#Concatenated simplified dynamics                  
#=============================================================
def Concatenated_simplified_dynamics(t0_fun, tf_fun, y0_fun, step_fun, Initial_dynamics_fun, tipping_points ,Dictionary_dynamics_fun ,hash_dynamics_fun,list_events,Terms_Ordered_var): 
    
    Simplified_dynamics_list=[];Simplified_times_list=[]

    Regime_results={}
    counter=0                                                                                       
    
    while t0_fun<tf_fun-step_fun:                                                                   
        print(Dictionary_dynamics_fun[hash_dynamics_fun])                                           
        
        #Solve individual dynamic
        #---------------------------------------------------------
        z_fun,time_fun,y_f_fun,type_event=Solve_Simplified_Model(y0_fun,t0_fun,tf_fun,tipping_points
,Dictionary_dynamics_fun[hash_dynamics_fun],[list_events[0], list_events[1]])                       
        #-----------------------------------------------------------

        print(" I am solving the dynamics")
        print("This is the solution of regime " + str(counter))
        print("Final bacteria")
        print(z_fun[0][-1])
        print("Final phage")
        print(z_fun[1][-1])
        print("Maximum bacteria")
        print(np.max(z_fun[0]))
        print("Maximum phage")
        print(np.max(z_fun[1]))
        print("Regime time")
        print(time_fun[-1] - time_fun[0])
        print("\n")

        Regime_results[counter]={'Bf':z_fun[0][-1],'Pf':z_fun[1][-1], 'MaxB':np.max(z_fun[0]), 'MaxP':np.max(z_fun[1]), 'tk_0':time_fun[0], 'tk_f': time_fun[-1]}

        #Update information for next iterations
        #-----------------------------------------------------------
        t0_fun=time_fun[-1];y0_fun=y_f_fun;counter+=1
        hash_dynamics_fun,Initial_dynamics_fun=Monitor_Dynamics(Initial_dynamics_fun,type_event,Terms_Ordered_var)
        #-----------------------------------------------------------
        
        #Store information of this iteration
        #--------------------------------------
        # if counter>1:
        #     z_funtest=[np.delete(z_fun[i], 0) for i in range(len(z_fun))]
        #     Simplified_dynamics_list.append(z_funtest)
        #     Simplified_times_list.append(time_fun[1:])
        # else:
        #     Simplified_dynamics_list.append(z_fun)
        #     Simplified_times_list.append(time_fun)

        Simplified_dynamics_list.append(z_fun)
        Simplified_times_list.append(time_fun)
            
        #--------------------------------------
        
        #------------------------------------------------------------
        
    Simplified_dynamics_fun=np.concatenate(Simplified_dynamics_list,axis=1)                   
    Simplified_times_fun=np.concatenate(Simplified_times_list)
    Critical_Times=[time[-1] for time in Simplified_times_list]
    #.............................................................
    
    return Simplified_dynamics_fun, Simplified_times_fun, Critical_Times, Regime_results
#=============================================================
    
#============================================================= 

#8. Find thresholds in a list
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Threshold_Finder(Vector, Threshold):
    
    Indices_Threshold=np.where(np.sign(Vector[:-1]-Threshold)!=\
    np.sign(Vector[1:] - Threshold) )[0] + 1
    
    return Indices_Threshold
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#9. Find Criticals
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Find_Criticals(df_fun, epsilon_value, dominant):
    
    #BURST
    #...............................................................
    #Critical indices for burst
    Per_Capita_Burst=df_fun['Per_Capita_Burst'].to_numpy()
    Critical_Ind_Burst=Threshold_Finder(Per_Capita_Burst, epsilon_value)
    #Critical concentrations
    CritConc_Burst=[];Critical_Times_Burst=[]
    for index in Critical_Ind_Burst:
        CritConc_Burst.append(df_fun.iloc[index]['Bacterial_Concentration'])
    #...............................................................

    #INFECTION
    #...............................................................
    #Critical indices 
    Per_Capita_Inf=df_fun['Per_Capita_Predation'].to_numpy()
    Critical_Ind_Inf=Threshold_Finder(Per_Capita_Inf, epsilon_value)
    #Critical concentration
    CritConc_Inf=[];Critical_Times_Infection=[]
    for index in Critical_Ind_Inf:
        CritConc_Inf.append(df_fun.iloc[index]['Phage_Concentration'])
    #...............................................................

    #SAVE DATA
    #...............................................................
    Critical_Indices={'Burst':Critical_Ind_Burst,'Predation':Critical_Ind_Inf}
    Critical_Concentrations={'Burst':CritConc_Burst,'Predation':CritConc_Inf}
    #...............................................................
    
    return Critical_Indices, Critical_Concentrations
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Find_one_per_ml(df_fun, column_bacteria_name, column_phage_name, Min_per_ml,threshold_2=0):
    
    #Bacteria
    #...............................................................
    #Critical indices for burst
    Bacteria=df_fun[column_bacteria_name].to_numpy()
    One_Bact_per_ml=Threshold_Finder(Bacteria, Min_per_ml)
    #Critical concentrations
    Bacteria_list=[];Critical_Times_Bacteria=[]
    for index in One_Bact_per_ml:
        Bacteria_list.append(df_fun.iloc[index][column_bacteria_name])
    #...............................................................

    #INFECTION
    #...............................................................
    #Critical indices 
    Phage=df_fun[column_phage_name].to_numpy()
    One_Phage_per_ml=Threshold_Finder(Phage, Min_per_ml)
    if threshold_2!=0:
        One_Phage_per_ml=Threshold_Finder(Phage, threshold_2)
    #Critical concentration
    Phage_list=[];Critical_Times_Phage=[]
    for index in One_Phage_per_ml:
        Phage_list.append(df_fun.iloc[index][column_phage_name])
    #...............................................................

    #SAVE DATA
    #...............................................................
    One_per_ml_Indices={'Bact':One_Bact_per_ml,'Phage':One_Phage_per_ml}
    One_per_ml_Concentrations={'Burst':Bacteria_list,'Phage':Phage_list}
    #...............................................................
    
    return One_per_ml_Indices, One_per_ml_Concentrations
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#10. Extract Active processes
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Get_Active_Bioprocesses(Dataframe,epsilon_fun):

    #New dataframe
    Binary_Dataframe=pd.DataFrame(columns=Dataframe.columns.tolist())
    #New dictionary
    Dict_Initial_Dynamics={}
    
    for process in Dataframe.columns.tolist():
        Binary_Dataframe[process]=(Dataframe[process] >=epsilon_fun).astype(int)
        Dict_Initial_Dynamics[process]=(Dataframe[process][0] >=epsilon_fun).astype(int)

    #Sum columns
    Binary_Dataframe['Total']=Binary_Dataframe.sum(axis=1)

    #Extract information
    Total_Active_Bioprocesses=Binary_Dataframe["Total"].values.tolist()
        
    return Total_Active_Bioprocesses, Dict_Initial_Dynamics
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#11. Logarithmic sampling of dataframe
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Logarithmic_Sampling_Vector(vector,Orders_of_Magnitude,First_Order):
    #Find where the vector crosses the first order of magnitude
    Threshold_First_Order=Threshold_Finder(vector,First_Order)
    Threshold_First_Order=Threshold_First_Order[0]

    #Fill the first order of sampling vector with all points from original vector
    First_Order_Vector=vector[:Threshold_First_Order]
    #All orders of magnitude have the same number of points
    Points_per_order=len(First_Order_Vector)

    #Run over higher order of magnitude and save log vectors
    Log_Sampling_Vectors=[]
    for Order in range(First_Order,Orders_of_Magnitude):
        #Find threshold of Nth order of magnitude
        Threshold_N_Order=Threshold_Finder(vector,10**Order)
        Threshold_N_Order=Threshold_N_Order[0]
        #Get original full vector between orders of magnitude
        N_Order_Vector=vector[Threshold_First_Order:Threshold_N_Order]
        #Sample original full vector with an increasing step
        Sampling_Step=int(0.1*(Points_per_order-1)*(10**(Order-1)))
        N_Order_Sampling=N_Order_Vector[0::Sampling_Step]
        #Save Nth order vector to list of lists 
        Log_Sampling_Vectors.append(N_Order_Sampling)

        #Define starting point for the next order of magnitude
        Threshold_First_Order=Threshold_N_Order

    #Sample vector with the remainings of the last order of magnitude
    #For instance, if your tf is 400 you dont have a full order of magnitude
    Last_Order=vector[Threshold_First_Order:]
    Last_Order_Sampling=Last_Order[0::Sampling_Step*10]

    #Concatenate all vectors of all orders of magnitude
    Log_Sampling_Vector_fun=np.concatenate((\
    First_Order_Vector,Log_Sampling_Vectors[0],Log_Sampling_Vectors[1],Last_Order_Sampling),axis=0)

    return Log_Sampling_Vector_fun
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#12 Total errors of simplified dynamics
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Get_Total_Error(True_Dynamics_Vec,Model_Dynamics_Vec, t_obs):

    #Error bacteria
    Bacteria_True=True_Dynamics_Vec[0];Bacteria_Model=Model_Dynamics_Vec[0]
    
    #Absolute errors
    Error_bact_vector=[]
    for (BT, BM) in zip(Bacteria_True, Bacteria_Model):
        Error_bact_i=(np.absolute(BT-BM))
        Error_bact_vector.append(Error_bact_i)

    #Error phage
    Phage_True=True_Dynamics_Vec[1];Phage_Model=Model_Dynamics_Vec[1]
    
    Error_phage_vector=[]
    for (PT, PM) in zip(Phage_True, Phage_Model):
        Error_phage_i=(np.absolute(PT-PM))
        Error_phage_vector.append(Error_phage_i)

    Total_error_bact=sum(Error_bact_vector)
    print(Total_error_bact)
    Total_error_phage=sum(Error_phage_vector)

    #Normalize by time
    #-------------------------------------
    Error_bact=Total_error_bact/t_obs
    Error_phage=Total_error_phage/t_obs
    #-------------------------------------
    
    return Error_bact, Error_phage
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#13. Relative errors of simplified dynamics
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Get_Relative_Error(True_Dynamics_Vec,Model_Dynamics_Vec):

    #Error bacteria
    Bacteria_True=True_Dynamics_Vec[0];Bacteria_Model=Model_Dynamics_Vec[0]

    Rel_Err_B=[]
    for (BT, BM) in zip(Bacteria_True, Bacteria_Model):
        Relative_Error_B_i=(np.absolute(BT-BM))/BT
        Rel_Err_B.append(Relative_Error_B_i)

    #Error phage
    Phage_True=True_Dynamics_Vec[1];Phage_Model=Model_Dynamics_Vec[1]
    
    Rel_Err_P=[]
    for (PT, PM) in zip(Phage_True, Phage_Model):
        Relative_Error_P_i=(np.absolute(PT-PM))/PT
        Rel_Err_P.append(Relative_Error_P_i)

    #Mean of both errors
    Rel_Error=[np.mean((bacteria,phage)) for (bacteria,phage) in zip(Rel_Err_B,Rel_Err_P)]
    
    #Mean total values
    #----------------------------------------------------------------
    Mean_Rel_Error_Bacteria=np.mean(Rel_Err_B)
    Mean_Rel_Error_Phage=np.mean(Rel_Err_P)
    Mean_Rel_Error=np.mean([Mean_Rel_Error_Phage,Mean_Rel_Error_Bacteria])
    #----------------------------------------------------------------
    
    return Rel_Error, Rel_Err_B, Rel_Err_P
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#14. Log-ratio error
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Get_logratio_Error(True_Dynamics_Vec,Model_Dynamics_Vec):

    #Log-ratio bacteria
    Bacteria_True=True_Dynamics_Vec[0];Bacteria_Model=Model_Dynamics_Vec[0]

    logratio_Err_B=[]
    for (BT, BM) in zip(Bacteria_True, Bacteria_Model):
        logratio_B_i=np.absolute(np.log(BM/BT))
        logratio_Err_B.append(logratio_B_i)

    #Log-ratio  phage
    Phage_True=True_Dynamics_Vec[1];Phage_Model=Model_Dynamics_Vec[1]
    
    logratio_Err_P=[]
    for (PT, PM) in zip(Phage_True, Phage_Model):
        logratio_P_i=np.absolute(np.log(PM/PT))
        logratio_Err_P.append(logratio_P_i)

    #Mean of both errors
    logratio_Error=[np.mean((bacteria,phage)) for (bacteria,phage) in zip(logratio_Err_B,logratio_Err_P)]
        
    return logratio_Error, logratio_Err_B, logratio_Err_P
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#11. Parameters, concentrations, and time for dominant timescales
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Initial_Configuration(Dominant_Timescale):
    #Scale independent parameters
    K=1e7;c=150;Volume=1
    if Dominant_Timescale=='r':
        r=0.9;m=2.8e-3;d=3e-8;
        y0=[1e3,1e4]
        Beq=m/(c*d);Peq=r/d
        time_0=0;time_f=14
        
    elif Dominant_Timescale=='m':
        r=3.7e-5;m=0.1;d=1e-10
        Beq=m/(c*d);Peq=r/d
        time_0=0.1;time_f=260
        y0=[10*Beq,1]
        
    elif Dominant_Timescale=='mr_Equ':
        d=3e-8;m=0.1
        r=m
        Beq=m/(c*d);Peq=r/d
        y0=[1.2*Beq,0.8*Peq]
        time_0=0;time_f=200

    elif Dominant_Timescale=='mr_Disequ':
        rmin=0.01;rmax=0.169;
        d=3e-8;m=0.1
        r=m
        Beq=m/(c*d);Peq=r/d
        y0=[2.5*Beq,Peq]
        time_0=0;time_f=200

    elif Dominant_Timescale=='r_Kick_Experiment':
        r=0.9;m=2.8e-3;d=3e-8;
        Beq=m/(c*d);Peq=r/d
        y0=[Beq,Peq]
        time_0=0;time_f=140

    elif Dominant_Timescale=='r_Kick_Experiment_Logistic':
        r=0.9;m=2.8e-3;d=3e-8;K=1e7
#        Bcrit=r/(c*d);Pcrit=r/d
        Beq=m/(c*d);Peq=(r/d)*(1-m/(K*c*d))
        y0=[Beq,Peq]
        time_0=0;time_f=int(2e5)

    elif Dominant_Timescale=='m_Kick_Experiment':
        r=3.7e-5;m=0.1;d=1e-10
        Beq=m/(c*d);Peq=(r/d)
        y0=[Beq,Peq]
        time_0=0;time_f=1000

    elif Dominant_Timescale=='m_Kick_Experiment_Logistic':
        r=3.7e-5;m=0.1;d=1e-10;K=1e7
        Beq=m/(c*d);Peq=(r/d)*(1-m/(K*c*d))
        y0=[Beq,Peq]
        time_0=0;time_f=int(2e5)
        
    elif Dominant_Timescale=='mr_Kick_Experiment':
        d=3e-8;m=0.1
        r=m
        Beq=m/(c*d);Peq=r/d
        y0=[Beq,Peq]
        time_0=0;time_f=1000

    elif Dominant_Timescale=='mr_Kick_Experiment_Logistic':
        d=3e-8;m=0.1;K=1e7
        r=m
        Beq=m/(c*d);Peq=r/d*(1-m/(K*c*d))
        y0=[Beq,Peq]
        time_0=0;time_f=int(2e5)
        
    Parameters={'r':r,'d':d,'c':c,'m':m,'K':K}

    return Parameters,y0,time_0,time_f
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##HYPERION FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#0.2. Split range of parameters in equispaced or logspaced intervals
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def LHS_intervals(min_value,max_value,Samples_fun,Logarithmic=0):
    interval=np.linspace(min_value, max_value, Samples_fun+1)
    
    if Logarithmic==1:
        interval=np.logspace(np.log10(min_value), np.log10(max_value),\
                             num=Samples_fun+1, endpoint=True, base=10)
        
    return interval 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.3. Latin_Hypercube_Sampling
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
#Ranges_Parameters: a dictionary where keys are names of parameters and values
#are a list of two elements: mimimum and maxium value of the interval
def LHS(Ranges_Parameters, Sampling_Points, Seed):
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

#0.4. Measurement Interval
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#This function gives you the experimental interval/span of time given a desired
#measurement time
def Measurement_Interval_fun(measurement_time,measure_error, time_vector):

    Initial_time=measurement_time - measure_error
    Final_time=measurement_time + measure_error

#Look for the chunk of time vector corresponding to the initial and final times
    Measure_Interval=[time_vector[np.abs(time_vector - Initial_time).argmin() :\
                                  np.abs(time_vector - Final_time).argmin()]][0]
    
    return Measure_Interval

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.5. Snapshots to Dataframe
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Given data of multiple measurement experiments, produce a dataframe with
#information about time experiments and measurements
def Snapshots_to_Dataframe(Experiments_Dict_var, Solutions_Dict_var, time_var, Measure_time_var):
    
    List_df=[]
    for Experiment_Number in Experiments_Dict_var:
        for Measurement in Experiments_Dict_var[Experiment_Number]:
            #Save snapshots of the model to list for Dataframe
            Times_Snapshot=Experiments_Dict_var[Experiment_Number][Measurement]
            counter_time=0
            for key1 in Solutions_Dict_var:
                for key2 in Solutions_Dict_var[key1]:
                    time_experiment=Times_Snapshot[counter_time]
                    index=np.where(time_var==time_experiment)
                    
                    List_df.append([Measure_time_var,Experiment_Number,Measurement
,time_var[index][0],Solutions_Dict_var[key1][key2][index][0],key2,key1])
                    counter_time+=1
                    
    Snapshot_df_var=pd.DataFrame(List_df,columns=["Theoretic Time","Experiment","Measurement","Time","Concentration","Bioagent","Type"])
    return Snapshot_df_var 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.6. Solutions to dataframe
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#This function converts a dictionary of solutions of a dynamic system of ODEs
#into a dataframe

def Solutions_to_Dataframe(Solutions_Dict_var,time_var):
    
    List_df=[]
    
    for key1 in Solutions_Dict_var:
        for key2 in Solutions_Dict_var[key1]:
            for snapshot in zip(Solutions_Dict_var[key1][key2],time_var):
                
                List_df.append([snapshot[1],snapshot[0],key2,key1])
                
    Solutions_df_var=pd.DataFrame(List_df,columns=[\
                "Time","Concentration", "Bioagent","Type"])
    
    return Solutions_df_var  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##OLD FUNCTIONS

#0.0. Solution to Lotka-Volterra
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Solve_Lotka_Volterra(Parameters, Initial_Conditions, time_vector, step):
    
    r=Parameters['r'];K=Parameters['K'];d=Parameters['d']
    c=Parameters['c'];m=Parameters['m']
    #0.1. Lotka-Volterra system
    #=================================================    
    def Lotka_Volterra(t,y):

        return [r*(1-(y[0]/float(K)))*y[0] - d*y[0]*y[1],
                c*d*y[0]*y[1] - m*y[1]]

    #Breaking condition for sensitive bacteria. Bottom
    def Min_Volume_B(t,y):
        return y[0] - 1
    Min_Volume_B.terminal=True
    Min_Volume_B.direction = 1
    #=================================================

    #0.2. MODEL SOLUTION
    #================================================================
    Events_Model=Min_Volume_B
    y0=[Initial_Conditions[0], Initial_Conditions[1] ]

    #0.2.1. Solver
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    sol_LV=solve_ivp(Lotka_Volterra,[time_vector[0],time_vector[-1]],y0,\
    method='RK45',dense_output=True,events=Events_Model,max_step=step)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.2.2. Change time vector in case event triggered
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    try:
        final_t=sol_LV.t_events[0][0]
        time=time_vector[:int(time_vector[-1])]
        print('Event found')
    except IndexError:
        pass
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    return sol_LV
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Solution to Marisa's experiment
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Solve_Experiment_Induction(Parameters, Initial_Conditions, time_vector, step):
    r=Parameters['r']
    K=Parameters['K']
    mu=Parameters['mu']
    c=Parameters['c']
    m=Parameters['m']

    #0.1. Induction system
    #=================================================
    def Experiment_Induction(t,y):
        
        return [r*(1-(y[0]/float(K)))*y[0] - mu*y[0],
                c*mu*y[0] - m*y[1]]
    
    #Breaking conditions for sensitive bacteria. Bottom
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    def Min_Volume_L(t,y):
        return y[0] - 1
    Min_Volume_L.terminal=False
    Min_Volume_L.direction = 1
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #=================================================

    #0.2. MODEL SOLUTION
    #================================================================
    Events_Model=[]
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
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    return sol_Exp
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
