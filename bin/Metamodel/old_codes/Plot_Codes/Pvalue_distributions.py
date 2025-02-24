#Script to plot pvalue distribution from the experiments of hyperion
import pandas as pd
import seaborn as sns
from Metamodel_Functions import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

#1. MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Input_Path='/home/sergio/work/Github/needle-finder/results/Data/Results_Marisa/'
Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/'


File_Name='Pvalues_36_Measurements_100_Number_of_Experiments.csv'
Pvalues_df=pd.read_csv(Input_Path + File_Name)


sns.set_context("talk") #paper, notebook, poster
g = sns.FacetGrid(Pvalues_df, col="Bioagent")
g.map_dataframe(sns.histplot, x="Pvalues", bins=33, log_scale=True).set(yscale='log')

g.axes[0,0].set_ylabel('Counts')
g.axes[0,1].set_xlabel('P-value')

plt.savefig(Output_Path+'Pvalue_Distribution.png',dpi=300)

plt.show()  








#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
