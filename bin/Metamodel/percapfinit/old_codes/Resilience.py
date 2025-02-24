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

                                                                                                  
from scipy.signal import find_peaks                                                                 
#++++++++++++++++++++++++++++++++++++++++++++




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
path='../../../data/input_output_data/'
file_name='Input_File_Kick.csv'

inputs=pd.read_csv(path+file_name)
inputs=inputs.set_index('description')                                                       
print(inputs)

scenarios_bacteria=['growth_decay_1e4_B_kick', 'growth_decay_5e3_B_kick', 'growth_decay_1e3_B_kick','growth_decay_5e2_B_kick', 'growth_decay_1e2_B_kick']

scenarios_phage=['growth_decay_1e6_P_kick', 'growth_decay_5e5_P_kick', 'growth_decay_1e5_P_kick','growth_decay_5e4_P_kick', 'growth_decay_1e4_P_kick']

scenarios_bacteria_phage=['growth_decay_1e4_B_1e6_P_kick', 'growth_decay_5e3_B_5e5_P_kick', 'growth_decay_1e3_B_1e5_P_kick', 'growth_decay_5e2_B_5e4_P_kick', 'growth_decay_1e2_B_1e4_P_kick']

all_scenarios=scenarios_bacteria+scenarios_phage+scenarios_bacteria_phage

Results=inputs.loc[all_scenarios]


Results['Delta_B_kick']=Results['B_eq']-Results['B0_kick']
Results['Delta_P_kick']=Results['P_eq']-Results['P0_kick']

#x axis
Results['Mean_Delta_B_kick_Delta_P_kick']=Results[['Delta_B_kick', 'Delta_P_kick']].mean(axis=1)
x=Results['Mean_Delta_B_kick_Delta_P_kick']

Results['Delta_B_min']=Results['B_eq-B_min']
Results['Delta_P_min']=Results['P_eq-P_min']

#y axis
Results['Mean_Delta_B_min_Delta_P_min']=Results[['Delta_B_min', 'Delta_P_min']].mean(axis=1)
y=Results['Mean_Delta_B_min_Delta_P_min']

print(Results[['B_eq','B0_kick','Delta_B_kick', 'P_eq','P0_kick','Delta_P_kick', 'Mean_Delta_B_kick_Delta_P_kick']])
print(x)
print(y)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

Output_Path='/home/sergio/work/Github/needle-finder/results/Plots/Paper_Matt/'

#Size of figures
cm = 1/2.54  # centimeters in inches                                                                
Width=8*cm;Height=5.0*cm #Width and height of plots                                                 

#Fontsizes
size_axis=7;size_ticks=6;size_title=5

#Extensions to save     
Extensions=['.svg', '.pdf']  

#Colors            
cmap='RdBu';cmap_pieces=matplotlib.cm.get_cmap(cmap)
Phage_color=cmap_pieces(0.1);Bacteria_color=cmap_pieces(0.9)                            

#Linewidth
width_line=2

#Figure
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
fig=figure(figsize=(Width, Height), dpi=300)

plt.scatter(x,y)
y_linear=x
plt.plot(x,y_linear,color='k',linestyle='dotted')


plt.xlabel('<B$_{eq}$ - B$_{k}$,P$_{eq}$ - P$_{k}>$',fontsize=size_axis)
plt.ylabel('<B$_{eq}$ - B$_{min}$, B$_{eq}$ - B$_{min}$> ',fontsize=size_axis)

# plt.xlabel('B$_{eq}$ - B$_{k}$',fontsize=size_axis)
# plt.ylabel('B$_{eq}$ - B$_{min}$ ',fontsize=size_axis)

plt.xticks(fontsize=size_ticks); plt.yticks(fontsize=size_ticks)  

sns.despine(top=True, right=True, left=False, bottom=False)

Name_Fig='Resilience'
[plt.savefig(Output_Path+Name_Fig+ext,dpi=300) for ext in Extensions]


plt.show()



