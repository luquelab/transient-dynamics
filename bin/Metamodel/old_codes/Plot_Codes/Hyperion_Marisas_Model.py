
#Import libraries
#+++++++++++++++++++++++++++++++++++++++
import json
import numpy as np
from scipy.integrate import odeint
import seaborn as sns
import math 
from scipy.integrate import solve_ivp
import pandas as pd
from decimal import Decimal 
from Metamodel_Functions import *
from matplotlib.pyplot import figure
import pickle
from scipy import stats
from statannotations.Annotator import Annotator
import matplotlib.ticker as mticker


pd.set_option('display.max_columns', None)

#++++++++++++++++++++++++++++++++

#0.Functions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Statistical Tests
#==============================================================================
def Statistical_Tests(Sample1, Sample2):
    Variance_1=Sample1.var(ddof=0)
    Variance_2=Sample2.var(ddof=0)
    
    if Variance_1!=Variance_2:
        Ttest=stats.ttest_ind(Sample1, Sample2, equal_var=False)

    else:
        Ttest=stats.ttest_ind(Sample1, Sample2)

    WMWtest=stats.mannwhitneyu(Sample1, Sample2)

    return Ttest, WMWtest
#==============================================================================


#0.2. Extract pseudo-experimental functions
#This function gives you measurements at experimental times each corresponding
#to a sampling from LHS
#=================================================================================
def LHS_Measurements_fun(Measurements_Model_var, Solutions_LHS_var):

    Total_Types=pd.unique(Measurements_Model_var["Total_Type"])
    Total_Type=Measurements_Model_var.iloc[0]["Total_Type"]
    
    #Extract Concentrations for all times and samplings from LHS
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #Filter by times
    Measurements_lhs_unfiltered=Solutions_LHS_var.loc[\
    Solutions_LHS_var['Time'].isin(Measurements_Model_var["Time"])]
    
    Measurements_lhs_unfiltered=Measurements_lhs_unfiltered[\
    Measurements_lhs_unfiltered["Total_Type"].isin(Total_Types)]
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Extract by sampling-time pairs. A sampling corresponds to a measurement
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #Define empty dataframe with same header
    Measurements_lhs_filtered=\
    pd.DataFrame(columns=Measurements_lhs_unfiltered.columns)

    for Pseudo_Experiment in zip(Measurements_Model_var["Time"],Measurements_Model_var["Measurement"],Measurements_Model_var["Total_Type"]):

        #Extract value that satisfies time and sampling conditions
        #................................................................
        Filter_by_time_and_sampling=Measurements_lhs_unfiltered[\
        (Measurements_lhs_unfiltered["Time"]==Pseudo_Experiment[0]) & \
        (Measurements_lhs_unfiltered["Sampling"]==Pseudo_Experiment[1]) &
        (Measurements_lhs_unfiltered["Total_Type"]==Pseudo_Experiment[2])]
        #................................................................

        #Append to dataframe
        #................................................................
        Measurements_lhs_filtered=\
        Measurements_lhs_filtered.append(Filter_by_time_and_sampling)
        #................................................................
        #-------------------------------------------------------------------------

    #Remove dummy index and reset
    #........................................................
    Measurements_lhs_filtered=\
    Measurements_lhs_filtered.drop(['Unnamed: 0'],axis=1)

    Measurements_lhs_filtered.reset_index(drop=True, inplace=True)
    #........................................................
    
    return Measurements_lhs_filtered
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#=================================================================================


#==============================================================================

#Boxplot
#==============================================================================
def Boxplot_fun(df_var, order_var, axes_var):

    #Plot box and swarmplot
    #---------------------------------------------------------------------------
    ax=sns.boxplot(x="Total_Type",y="Concentration",data=df_var,order=order_var)
    ax=sns.swarmplot(x="Total_Type",y="Concentration", data=df_var,color=".25")
    #---------------------------------------------------------------------------

    #Labels
    #---------------------------------------------------------------------------
    ax.set_xlabel('',fontsize=axes_var['font']['axis'])
    ax.set_ylabel(axes_var['labels']['y'],fontsize=axes_var['font']['axis'])
    #---------------------------------------------------------------------------
    #Ticks
    #---------------------------------------------------------------------------
    ax.set_xticklabels(['Wild type', 'Mutant'])
    ax.tick_params(axis='both', which='major',labelsize=axes_var['font']['ticks'])

    #Scientific notation for ticks
    f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
    plt.gca().yaxis.set_major_formatter(mticker.FuncFormatter(g))
    #---------------------------------------------------------------------------
    
    return ax
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Load data
#=========================================================================
#Pickles
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Data_Path='/home/sergio/work/Github/needle-finder/results/python_pickles/'
Data=pickle.load(open(Data_Path + 'Pickle_Test.p', 'rb') )
Solution_wt=Data[0];Solution_sp=Data[1];time=Data[2]
Measure_Err=Data[4];
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/'

#Dataframe of experiments
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Data_Path_1='/home/sergio/work/Github/needle-finder/results/Data/Results_Marisa/'
Measurement_Time=10
Snapshot_Experiments_df=pd.read_csv(Data_Path_1 + 'Snapshot_Experiment_time' + str(Measurement_Time) + '.csv')
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Latin Hypercube Sampling
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Interval_Width=0.1
Solutions_LHS=pd.read_csv(Data_Path_1 + 'Solutions_LHS_'+\
            str(Interval_Width)+'_Interval_variation.csv')

Solutions_LHS["Total_Type"] = \
Solutions_LHS["Bioagent"] + Solutions_LHS["Type"]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#=========================================================================


#Take single experiment without and with LHS
#=========================================================================
Single_Experiment_df=Snapshot_Experiments_df[\
Snapshot_Experiments_df['Theoretic Time']==10.0]
Single_Experiment_df=Single_Experiment_df.drop(['Experiment'],axis=1)
Single_Experiment_df=Single_Experiment_df.drop(['Unnamed: 0'],axis=1)
#=========================================================================

#Extract VMR
#=========================================================================
Single_Experiment_df["Total_Type"] = \
Single_Experiment_df["Bioagent"] + Single_Experiment_df["Type"]
    
Cells_df=Single_Experiment_df[Single_Experiment_df["Bioagent"]=='CFU']
Cells_df.reset_index(drop=True, inplace=True)

Cells_df_lhs=LHS_Measurements_fun(Cells_df,Solutions_LHS)

Phages_df=Single_Experiment_df[Single_Experiment_df["Bioagent"]=='PFU']
Phages_df.reset_index(drop=True, inplace=True)

Phages_df_lhs=LHS_Measurements_fun(Phages_df,Solutions_LHS)


VMR=Phages_df['Concentration']/Cells_df['Concentration']
VMR_series = pd.Series(VMR, name="Concentration")
VMR_df=VMR_series.to_frame()
VMR_df['Bioagent']='VMR'
VMR_df['Measurement']=Cells_df['Measurement']
VMR_df['Time']=Cells_df['Time']
VMR_df['Type']=Cells_df['Type']
VMR_df['Total_Type']=VMR_df['Bioagent'] + VMR_df['Type']

VMR_lhs=Phages_df_lhs['Concentration']/Cells_df_lhs['Concentration']
VMR_lhs_series = pd.Series(VMR_lhs, name="Concentration")
VMR_lhs_df=VMR_lhs_series.to_frame()
VMR_lhs_df['Bioagent']='VMR'
VMR_lhs_df['Sampling']=Cells_df_lhs['Sampling']
VMR_lhs_df['Time']=0.5*(Cells_df_lhs['Time'] + Phages_df_lhs['Time'])
VMR_lhs_df['Type']=Cells_df_lhs['Type']
VMR_lhs_df['Total_Type']=VMR_lhs_df['Bioagent'] + VMR_lhs_df['Type']
#=========================================================================

#Statistical Tests
#=========================================================================

#Filter dataframes by biological agent
#---------------------------------------------------------------------------------

#Extract pseudo experimental data from CFU wt
#.................................................................................
#Extract data from the solution of the model
Measurements_CFU_wt=Single_Experiment_df[Single_Experiment_df['Total_Type']=='CFUwt']
#Change model data by LHS replicates
Measurements_CFU_wt_lhs=Cells_df_lhs[Cells_df_lhs['Total_Type']=='CFUwt']
#.................................................................................

#Extract pseudo experiment data from CFU sp
#.................................................................................
Measurements_CFU_sp=Single_Experiment_df[Single_Experiment_df['Total_Type']=='CFUsp']
#Change model data by LHS replicates
Measurements_CFU_sp_lhs=Cells_df_lhs[Cells_df_lhs['Total_Type']=='CFUsp']
#.................................................................................

#Extract pseudo experiment data from PFU wt
#.................................................................................
Measurements_PFU_wt=Single_Experiment_df[Single_Experiment_df['Total_Type']=='PFUwt']
#Change model data by LHS replicates
Measurements_PFU_wt_lhs=Phages_df_lhs[Phages_df_lhs['Total_Type']=='PFUwt']
#.................................................................................

#Extract pseudo experiment data from PFU sp
#.................................................................................
Measurements_PFU_sp=Single_Experiment_df[Single_Experiment_df['Total_Type']=='PFUsp']
#Change model data by LHS replicates
Measurements_PFU_sp_lhs=Phages_df_lhs[Phages_df_lhs['Total_Type']=='PFUsp']
#................................................................................
#---------------------------------------------------------------------------------

#Calculate VMR measurements from CFU and PFU measurements
#---------------------------------------------------------------------------------
#Wild Type
Measurements_VMR_wt_lhs=VMR_lhs_df[VMR_lhs_df['Total_Type']=='VMRwt']
#sp
Measurements_VMR_sp_lhs=VMR_lhs_df[VMR_lhs_df['Total_Type']=='VMRsp']
#---------------------------------------------------------------------------------

#.................................................................................

#---------------------------------------------------------------------------------

#Statistical tests
#--------------------------------------------------------------------------
#Microbial
Ttest_Cells_lhs, WMWtest_Cells_lhs=Statistical_Tests(Measurements_CFU_wt_lhs['Concentration'],Measurements_CFU_sp_lhs['Concentration'])

Ttest_Cells, WMWtest_Cells=Statistical_Tests(Measurements_CFU_wt['Concentration'],Measurements_CFU_sp['Concentration'])

#Phage
Ttest_Phages_lhs,WMWtest_Phages_lhs=Statistical_Tests(Measurements_PFU_wt_lhs['Concentration'], Measurements_PFU_sp_lhs['Concentration'])

Ttest_Phages,WMWtest_Phages=Statistical_Tests(Measurements_PFU_wt['Concentration'], Measurements_PFU_sp['Concentration'])

#VMR
VMR_wt=VMR_df[VMR_df['Total_Type']=='VMRwt']
VMR_sp=VMR_df[VMR_df['Total_Type']=='VMRsp']

Ttest_VMRs_lhs, WMWtest_VMRs_lhs=Statistical_Tests(Measurements_VMR_wt_lhs['Concentration'], Measurements_VMR_sp_lhs['Concentration'])

Ttest_VMRs,WMWtest_VMRs=Statistical_Tests(VMR_wt['Concentration'],VMR_sp['Concentration'])
#================================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#Plots
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Parameters and fontsizes
#========================
size_axis=18
size_ticks=16
size_letter=15
#========================

#Plot 1 - Bacteria and phages
#===============================================================================
#------------------------------------------------------------------------------
plt.plot(time, Solution_wt.sol(time).T, linewidth=3)
plt.plot(time, Solution_sp.sol(time).T, linewidth=3)
plt.axvline(x=Measurement_Time, color='k')

plt.axvspan(Measurement_Time-Measure_Err, Measurement_Time+Measure_Err,alpha=0.5, color='gray')

sns.scatterplot(data=Single_Experiment_df, x="Time", y="Concentration",s=150,color='k')
#------------------------------------------------------------------------------

#Axes limits
#--------------------------------------------
Max_wt=np.max( Solution_wt.sol(time).T )
Max_sp=np.max( Solution_sp.sol(time).T )
Max_conc=np.max([Max_wt, Max_sp])
plt.xlim([time[0], time[-1]])
#plt.xlim([9.25, 10.75])
plt.ylim([1,1.25*Max_conc ])
#--------------------------------------------

#Labels, legend and scale
#-------------------------------------------------------------------------
plt.xlabel(r'time $(h)$',fontsize=size_axis)
plt.ylabel(r'Concentration $(ml^{-1})$',fontsize=size_axis)
plt.legend([r'B_WT', 'Ph_WT', r'B_spoT-', 'Ph_SpoT-' ],loc='upper right')
plt.yscale('log')
#-------------------------------------------------------------------------

plt.show()
#=========================================================================


#Plot 2 - VMR
#===============================================================================
#------------------------------------------------------------------------------
VMR_wt_solution=np.divide(Solution_wt.sol(time)[1].T,Solution_wt.sol(time)[0].T)
VMR_sp_solution=np.divide(Solution_sp.sol(time)[1].T,Solution_sp.sol(time)[0].T)

plt.plot(time, VMR_wt_solution, color='blue',linewidth=3)
plt.plot(time, VMR_sp_solution, color='green',  linewidth=3)

plt.axvline(x=Measurement_Time, color='k')
plt.axvspan(Measurement_Time-Measure_Err, Measurement_Time+Measure_Err,alpha=0.5, color='gray')

sns.scatterplot(data=VMR_df, x="Time", y="Concentration", style="Measurement",s=100,color='k')
#------------------------------------------------------------------------------

#Labels, legend and scale
#-------------------------------------------------------------------------
plt.xlabel(r'time $(h)$',fontsize=size_axis)
plt.ylabel(r'VMR $(Phages/Cells)$',fontsize=size_axis)
plt.legend([r'VMR_wt', 'VMR_sp'],loc='upper right')
plt.yscale('log')
#-------------------------------------------------------------------------

plt.show()
#=========================================================================

#--------------------------------------------------------------------------

#Bacterial boxplot
#---------------------------------------------------------------------------------
axes_cells={\
'labels':{'x':'', 'y':'Concentrations (cells/ml)'},\
'ticks': {'x':['Wild_type', 'Mutant']},
'font':  {'axis':size_axis, 'ticks':size_ticks}}

figure(figsize=(12, 8), dpi=100)

order=["CFUwt", "CFUsp"]
Boxplot_cells=Boxplot_fun(Cells_df, order, axes_cells)

Boxplot_cells.set_ylim([5e8, 1.1e9 ])

#Annotate results
#.................................................................................
Boxplot_cells.text(0.75, 8e8, "WMW: "+ str(WMWtest_Cells[0]) + '\n' +\
"p-value: "+ str(WMWtest_Cells[1])+ '\n' + '\n' +\
"T-test: " + str(Ttest_Cells[0])+'\n'+"p-value: "+str(Ttest_Cells[1]),\
bbox=dict(facecolor='red', alpha=0.5),fontsize=size_ticks)
#.................................................................................

plt.savefig(Output_Path+'Boxplot_Cells.png',dpi=300) 

plt.show()
#---------------------------------------------------------------------------------

#Phage boxplot
#---------------------------------------------------------------------------------
axes_phage={\
'labels':{'x':'', 'y':'Concentrations (phages/ml)'},\
'ticks': {'x':['Wild_type', 'Mutant']},
'font':  {'axis':size_axis, 'ticks':size_ticks}}

figure(figsize=(12, 8), dpi=100)

order=["PFUwt", "PFUsp"]
Boxplot_phages=Boxplot_fun(Phages_df, order, axes_phage)

#Annotate
#.................................................................................
Boxplot_phages.text(0.75, 11, "WMW: "+ str(WMWtest_Phages[0]) + '\n' +\
"p-value: "+ str(WMWtest_Phages[1])+ '\n' + '\n' +\
"T-test: " + str(Ttest_Phages[0])+'\n'+"p-value: "+str(Ttest_Phages[1]),\
bbox=dict(facecolor='red', alpha=0.5),fontsize=size_ticks)
#.................................................................................

plt.savefig(Output_Path+'Boxplot_phages.png',dpi=300) 

plt.show()
#---------------------------------------------------------------------------------

#VMR boxplot
#---------------------------------------------------------------------------------
axes_VMR={\
'labels':{'x':'', 'y':'VMR (phages/cell)'},\
'ticks': {'x':['Wild_type', 'Mutant']},
'font':  {'axis':size_axis, 'ticks':size_ticks}}

figure(figsize=(12, 8), dpi=100)

order=["VMRwt", "VMRsp"]
Boxplot_VMR=Boxplot_fun(VMR_df, order, axes_VMR)

#Annotate
#.................................................................................
Boxplot_VMR.text(0.75, 1e-8, "WMW: "+ str(WMWtest_VMRs[0]) + '\n' +\
"p-value: "+ str(WMWtest_VMRs[1])+ '\n' + '\n' +\
"T-test: " + str(Ttest_VMRs[0])+'\n'+"p-value: "+str(Ttest_VMRs[1]),\
bbox=dict(facecolor='red', alpha=0.5),fontsize=size_ticks)
#.................................................................................

Boxplot_VMR.set_ylim([0.2e-8, 1.9e-8 ])

plt.savefig(Output_Path+'Boxplot_VMRs.png',dpi=300) 
plt.show()
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

#Scatterplot Cells VLPs
#---------------------------------------------------------------------------------
#fig=figure(figsize=(8, 6), dpi=80)  
fig,ax_vlp = plt.subplots()

ax_vlp.scatter(Measurements_CFU_wt['Concentration'], Measurements_PFU_wt['Concentration'],color='red',marker='^',label='wild_type_vlp')

ax_vlp.scatter(Measurements_CFU_sp['Concentration'], Measurements_PFU_sp['Concentration'],color='red', marker='s',label='mutant_vlp')

ax_vlp.set_yscale('log')

ax_vmr=ax_vlp.twinx()

ax_vmr.scatter(Measurements_CFU_wt['Concentration'], VMR_wt['Concentration'],color='blue', marker = '^',label='wild_type_vmr')
ax_vmr.scatter(Measurements_CFU_sp['Concentration'], VMR_sp['Concentration'],color='blue', marker='s',label='mutant_vmr')


ax_vlp.set_xlabel('Cells/ml',fontsize=size_axis);
ax_vlp.set_ylabel('phages/ml',fontsize=size_axis)
ax_vmr.set_ylabel('VMR',fontsize=size_axis)

ax_vmr.set_yscale('log')

fig.legend(loc="upper left")

plt.show()
#---------------------------------------------------------------------------------
