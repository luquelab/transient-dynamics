
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




#Load data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Data_Path_1='/home/sergio/work/Github/needle-finder/results/Data/Results_Marisa/'
Interval_Width=0.1
LHS_Solutions_df=pd.read_csv(Data_Path_1 + 'Solutions_LHS_'+\
                str(Interval_Width)+'_Interval_variation.csv')

Solutions_df=pd.read_csv(Data_Path_1 + 'Solutions.csv')
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Plot
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Plots
#=========================================================================
LHS_Solutions_df["Total_Type"]=LHS_Solutions_df["Bioagent"] + LHS_Solutions_df["Type"]
Solutions_df["Total_Type"]=Solutions_df["Bioagent"] + Solutions_df["Type"]
ax_LHS=sns.lineplot(data=LHS_Solutions_df, x="Time", y="Concentration", hue="Total_Type")
ax=sns.lineplot(data=Solutions_df, x="Time", y="Concentration", hue="Total_Type")
ax.lines[0].set_linestyle("--")
ax.lines[1].set_linestyle("--")
ax.lines[2].set_linestyle("--")
ax.lines[3].set_linestyle("--")
#=========================================================================

#Parameters and fontsizes
#========================
size_axis=18
size_ticks=16
size_letter=15
#========================

#Logscale
#========================
ax_LHS.set_yscale('log')
ax.set_yscale('log')
#========================

plt.show()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
