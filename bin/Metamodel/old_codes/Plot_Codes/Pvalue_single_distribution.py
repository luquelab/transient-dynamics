
#Import libraries
#+++++++++++++++++++++++++++++++++++++++
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
from statannot import add_stat_annotation
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


#Boxplot
#==============================================================================
def Boxplot(dataframe, Statistics):

    ax=sns.boxplot(x="Total_Type",y="Concentration", data=dataframe)
    ax=sns.swarmplot(x="Total_Type",y="Concentration", data=dataframe,color=".25")

    add_stat_annotation(ax,data=dataframe, x="Total_Type", y="Concentration", box_pairs=[("CFUwt", "CFUsp")], test='Mann-Whitney', text_format='star', loc='outside', verbose=2)

    ax.text(0.75, 8e8, "U1: "+ str(WMWtest_Cells[0]) + '\n' +\
    "p-value: "+ str(WMWtest_Cells[1]) , bbox=dict(facecolor='red', alpha=0.5))

    ax.set_xlabel('',fontsize=size_axis);
    ax.set_ylabel('Concentrations (cells/ml)',fontsize=size_axis)

    ax.set_xticklabels(['Wild type', 'Mutant'])
    ax.tick_params(axis='both', which='major', labelsize=size_ticks)
    
    return None
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Load data
#=========================================================================
Data_Path='/home/sergio/work/Github/needle-finder/results/python_pickles/'
Data=pickle.load(open(Data_Path + 'Pickle_Test.p', 'rb') )
Solution_wt=Data[0];Solution_sp=Data[1];time=Data[2]
Measure_Time=Data[3];Measure_Err=Data[4];Measure_Experiments_plot_df=Data[5];
#=========================================================================

#Extract VMR
#=========================================================================
Measure_Experiments_plot_df["Total_Type"] = \
Measure_Experiments_plot_df["Bioagent"] + Measure_Experiments_plot_df["Type"]

Cells_df=Measure_Experiments_plot_df[\
Measure_Experiments_plot_df["Bioagent"]=='CFU']
Cells_df.reset_index(drop=True, inplace=True)

Phages_df=Measure_Experiments_plot_df[\
Measure_Experiments_plot_df["Bioagent"]=='PFU']
Phages_df.reset_index(drop=True, inplace=True)

VMR=Phages_df['Concentration']/Cells_df['Concentration']
VMR_series = pd.Series(VMR, name="Concentration")
VMR_df=VMR_series.to_frame()
VMR_df['Experiment']=Cells_df['Experiment']
VMR_df['Measurement']=Cells_df['Measurement']
VMR_df['Time']=Cells_df['Time']
VMR_df['Bioagent']='VMR'
VMR_df['Type']=Cells_df['Type']
VMR_df['Total_Type']=VMR_df['Bioagent'] + VMR_df['Type']
#=========================================================================

#Statistical Tests
#=========================================================================

#Filter dataframes by biological agent
#--------------------------------------------------------------------------

Measurements_CFU_wt=Measure_Experiments_plot_df[Measure_Experiments_plot_df['Total_Type']=='CFUwt']
Measurements_CFU_sp=Measure_Experiments_plot_df[Measure_Experiments_plot_df['Total_Type']=='CFUsp']
Measurements_PFU_wt=Measure_Experiments_plot_df[Measure_Experiments_plot_df['Total_Type']=='PFUwt']
Measurements_PFU_sp=Measure_Experiments_plot_df[Measure_Experiments_plot_df['Total_Type']=='PFUsp']

VMR_wt=VMR_df[VMR_df['Total_Type']=='VMRwt']
VMR_sp=VMR_df[VMR_df['Total_Type']=='VMRsp']


Number_of_Experiments=Measure_Experiments_plot_df.Experiment.unique()
print(Number_of_Experiments)
Pvalues_list=[]

for Experiment in Number_of_Experiments:

    #-----------------------------------------------------------------------------
    Measurements_CFU_wt_Experiment=Measurements_CFU_wt[Measurements_CFU_wt['Experiment']==Experiment]

    Measurements_CFU_sp_Experiment=Measurements_CFU_sp[Measurements_CFU_sp['Experiment']==Experiment]

    Ttest_Cells,WMWtest_Cells=Statistical_Tests(\
    Measurements_CFU_wt_Experiment['Concentration'],Measurements_CFU_sp_Experiment['Concentration'])
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    Measurements_PFU_wt_Experiment=Measurements_PFU_wt[Measurements_PFU_wt['Experiment']==Experiment]
    
    Measurements_PFU_sp_Experiment=Measurements_PFU_sp[Measurements_PFU_sp['Experiment']==Experiment]

    Ttest_Phages,WMWtest_Phages=Statistical_Tests(\
    Measurements_PFU_wt_Experiment['Concentration'],Measurements_PFU_sp_Experiment['Concentration'])
    #-----------------------------------------------------------------------------
    
    #-----------------------------------------------------------------------------
    VMR_wt_Experiment=VMR_wt[VMR_wt['Experiment']==Experiment]
    VMR_sp_Experiment=VMR_sp[VMR_sp['Experiment']==Experiment]
    Ttest_VMRs, WMWtest_VMRs=Statistical_Tests(VMR_wt_Experiment['Concentration'], VMR_sp_Experiment['Concentration'])
    #-----------------------------------------------------------------------------
    Pvalues_list.append([WMWtest_Cells[1],WMWtest_Phages[1],WMWtest_VMRs[1]])


    
Pvalues_df=pd.DataFrame(Pvalues_list, columns=["Cells","Phages","VMRs"])
print(Pvalues_df)

#Parameters and fontsizes
#======================== 
size_axis=18
size_ticks=16
size_letter=15
#======================== 



ax=sns.histplot(data=Pvalues_df, x="VMRs", bins=10)

ax.set_xlabel('VMRs',fontsize=size_axis)
plt.show()
#--------------------------------------------------------------------------
