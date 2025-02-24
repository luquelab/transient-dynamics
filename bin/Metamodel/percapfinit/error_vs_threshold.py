#25_03_2024. This code represents the error vs the threshold of active mechanisms.

import math
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib                                                    
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
from Metamodel_Functions import *

#===========================================
def percentage_fun(list_error):
    list_error=[100*i for i in list_error]

    return list_error
#===========================================

#Read dataset
path='../../../data/input_output_data/' 
file_name='errors_new.csv' 
data=pd.read_csv(path+file_name)
print(data['category'])

print(data[ (data['threshold']==1.0) & (data['category']=='max' )] )

print("\n")
print(data[ (data['category']=='median') | (data['category']=='mad') ] )
median_errors=data[ (data['category']=='median') | (data['category']=='mad') ]
median_errors.to_csv('../../../data/median_errors.csv', index=True)

#scenarios to plot: growth,decay,growth-decay,growth-decay_regimes,no-growth-no-decay
scenario_to_plot='growth'

#Error type
err_type='logratio' #logratio, relative

data=data[data['scenario']==scenario_to_plot]
data=data.sort_values(by=['threshold'])

max_errors=data[data['category']=='max']
end_errors=data[data['category']=='end']

median_errors=data[data['category']=='median']

median_absolute_deviation=data[data['category']=='mad']

mean_errors=data[data['category']=='mean']
std_mean_errors=data[data['category']=='std']
print(mean_errors)
print(std_mean_errors)

#Figure settings
#===================================================================
#Path to save figure                                                
Output_Path='../../../results/Plots/plots_paper/'                   

#Fontsizes                                                          
size_axis=7;size_ticks=6;size_title=5                               

#Figure Size                                                        
cm = 1/2.54  # centimeters in inches                                
Width=8*cm;Height=4*cm #Width and height of plots                   

#Gridspec parameters                                                
Rows=1;Cols=2                                                       

#Colors                                                             
cmap='RdBu';cmap_pieces= plt.get_cmap(cmap)               
Phage_color=cmap_pieces(0.1);Bacteria_color=cmap_pieces(0.9)
Phage_color_end=cmap_pieces(0.3);Bacteria_color_end=cmap_pieces(0.7)

#Extensions to save                                                 
Extensions=['.svg','.pdf']                                          
#Linewidth                                                          
width_line=1                                                        
#================================================================

thresholds=max_errors['threshold'].tolist()

#Max and end figures
#--------------------------------------------------------------------
fig=figure(figsize=(Width, Height), dpi=300)

bacteria_max=max_errors['bacteria' + '_' + err_type]
#bacteria_max=percentage_fun(bacteria_max)
phage_max=max_errors['phage' + '_' + err_type]
#phage_max=percentage_fun(phage_max)

bacteria_end=end_errors['bacteria' + '_' + err_type]
#bacteria_end=percentage_fun(bacteria_end)
phage_end=end_errors['phage' + '_' + err_type]
#phage_end=percentage_fun(phage_end)

#Plot data
plt.plot(thresholds,bacteria_max,'-.^',markersize=3,color=Bacteria_color,linewidth=1,label='Prey max')
plt.plot(thresholds,phage_max,'-.^',markersize=3,color=Phage_color,linewidth=1,label='Predator max')

plt.plot(thresholds,bacteria_end,':*',markersize=3,color=Bacteria_color_end,linewidth=1,label='Prey end')
plt.plot(thresholds,phage_end,':*',markersize=3,color=Phage_color_end,linewidth=1,label='Predator end')

#Ticks, scale, limits and legend
xticks_labels=[0.1, 0.5, 1, 1.5, 2]
plt.xticks(xticks_labels, fontsize=size_ticks)
#yticks_labels=[10$^, 0.5, 1, 1.5, 2]
plt.yticks(fontsize=size_ticks)

plt.yscale('log');
plt.ylim([1e-4,10])


#Capitalize error name and add string
ylabel=err_type[0].upper() + err_type[1:] + '  error'

plt.ylabel(ylabel,fontsize=size_axis)
plt.xlabel('Critical threshold',fontsize=size_axis)

plt.legend(loc='best',fontsize=size_ticks,frameon=False)

Name_Figure='critical_threshold_end_max_'+err_type+'_'+str(scenario_to_plot)
for ext in Extensions:                                               
    plt.savefig(Output_Path+Name_Figure+ext,dpi=300) 

plt.show()
#--------------------------------------------------------------------

#Median
#--------------------------------------------------------------------
fig=figure(figsize=(Width, Height), dpi=300)

median_bacteria=median_errors['bacteria' + '_' + err_type]
#median_bacteria=percentage_fun(median_bacteria)
median_phage=median_errors['phage' + '_' + err_type]
#median_phage=percentage_fun(median_phage)

mad_bacteria=median_absolute_deviation['bacteria' + '_' + err_type]
#mad_bacteria=percentage_fun(mad_bacteria)
mad_phage=median_absolute_deviation['phage' + '_' + err_type]
#mad_phage=percentage_fun(mad_phage)



plt.errorbar(thresholds,median_bacteria, yerr=mad_bacteria,color=Bacteria_color,linewidth=1,marker='o',markersize=3, linestyle=':', label='Prey median')

plt.errorbar(thresholds,median_phage, yerr=mad_phage,color=Phage_color,linewidth=1,marker='o',markersize=3,linestyle=':',label='Predator median')

#Ticks, scale, limits and legend
xticks_labels=[0.1, 0.5, 1, 1.5, 2]
plt.xticks(xticks_labels, fontsize=size_ticks)
plt.yticks(fontsize=size_ticks)

plt.yscale('log');
plt.ylim([1e-4,10])

 

#plt.ylabel(err_type + ' (%)',fontsize=size_axis)
plt.ylabel(ylabel,fontsize=size_axis)
plt.xlabel('Critical threshold',fontsize=size_axis)

plt.legend(loc='best',fontsize=size_ticks,frameon=False)

Name_Figure='critical_threshold_median_'+err_type+'_'+str(scenario_to_plot)
for ext in Extensions:                                               
    plt.savefig(Output_Path+Name_Figure+ext,dpi=300) 

plt.show()
#--------------------------------------------------------------------


#Mean and std figure
#----------------------------------------------------------------
mean_bacteria=mean_errors['bacteria' + '_' + err_type]
mean_phage=mean_errors['phage' + '_' + err_type]
#mean_bacteria=percentage_fun(mean_bacteria)
#mean_phage=percentage_fun(mean_phage)

std_bacteria=std_mean_errors['bacteria' + '_' + err_type]
#std_bacteria=percentage_fun(std_bacteria)
std_phage=std_mean_errors['phage' + '_' + err_type]
#std_phage=percentage_fun(std_phage)


fig=figure(figsize=(Width, Height), dpi=300)

plt.errorbar(thresholds,mean_bacteria, yerr=std_bacteria,color=Bacteria_color,linewidth=1,marker='o',markersize=3,label='Prey mean')

plt.errorbar(thresholds,mean_phage, yerr=std_phage,color=Phage_color,linewidth=1,marker='o',markersize=3,label='Predator mean')

#Ticks, scale, limits and legend
xticks_labels=[0.1, 0.5, 1, 1.5, 2]
plt.xticks(xticks_labels, fontsize=size_ticks)
plt.yticks(fontsize=size_ticks)

plt.yscale('log');
plt.ylim([1e-4,10])
#plt.ylim([1e-2,200])

#plt.ylabel(err_type + ' (%)',fontsize=size_axis)
plt.ylabel(ylabel ,fontsize=size_axis)
plt.xlabel('Critical threshold',fontsize=size_axis)

plt.legend(loc='best',fontsize=size_ticks,frameon=False)

Name_Figure='critical_threshold_mean_'+err_type+'_'+str(scenario_to_plot)
for ext in Extensions:                                               
    plt.savefig(Output_Path+Name_Figure+ext,dpi=300) 

plt.show()
#-----------------------------------------------------------------

