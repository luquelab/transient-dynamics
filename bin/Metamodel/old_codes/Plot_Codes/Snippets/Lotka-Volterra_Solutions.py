#2021-10-22. A snippet for a simple lotka-volterra plot.
#2022-06-08. It includes an event in the solution, vertical line plot and xtick positions and labels to illustrate the differences

#Libraries
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy.integrate import solve_ivp

#Lotka-Volterra equation
def Lotka_Volterra(t,y):
    return [r*y[0] - d*y[0]*y[1],
           c*d*y[0]*y[1] - m*y[1]]
def Extinction_prey(t,y):
    return y[0] -1
Extinction_prey.terminal=False
Extinction_prey.direction=0

#Parameters,time and solver
#====================================================================
r=0.9;m=2.8e-3;d=3e-8;c=150;
y0=[1e3,1e4]
time_0=0;time_f=15
Resolution=100;Steps=time_f*Resolution + 1
time=np.linspace(time_0,time_f,Steps)
step=0.01
raw_solution_LV=solve_ivp(Lotka_Volterra,[time[0],time[-1]],y0,\
method='RK45',dense_output=True,events=Extinction_prey,max_step=step)
extinction_time=raw_solution_LV.t_events[0][0]
plot_solution_LV=raw_solution_LV.sol(time)
prey=plot_solution_LV[0].T;predator=plot_solution_LV[1].T
#====================================================================

#Path to save figure
Parent_Path='/home/sergio/work/Github/needle-finder/'
Output_Path=Parent_Path+'results/Plots/'

#Fontsizes
size_axis=7;size_ticks=5;size_title=5 

#Figure Size
cm = 1/2.54  # centimeters in inches
Width=8*cm;Height=4*cm #Width and height of plots
fig=figure(figsize=(Width, Height), dpi=300)

#Colors
cmap='RdBu';cmap_pieces= matplotlib.cm.get_cmap(cmap)
predator_col=cmap_pieces(0.1);prey_col=cmap_pieces(0.9) 

#Plot Figure
plt.plot(time,prey,color=prey_col,linestyle='solid',linewidth=1,label='Prey')
plt.plot(time,predator,color=predator_col,linestyle='solid',linewidth=1,label='Predator')

#Vertical line
plt.axvline(extinction_time,color='k',linewidth=1,linestyle='dotted')

#Axes, title and ticks
plt.ylabel('Concentration (cells/ml)',fontsize=size_axis)
plt.xlabel('time (h)',fontsize=size_axis)
title_variable='Lotka-Volterra'
plt.title('Title of %s figure'  %title_variable,fontsize=size_title)
xticks_pos=[0,4,8,extinction_time,14]
xticks_labels=[0,4,8,'Extinction',14]
plt.xticks(xticks_pos,xticks_labels,fontsize=size_ticks);plt.yticks(fontsize=size_ticks)


#axes limits
plt.xlim(xmin=time_0,xmax=time_f)
ymin=min(min(prey),min(predator))
ymax=max(max(prey),max(predator))
                            
plt.ylim(ymin=1,ymax = 100*ymax)

#Legend and scale
plt.legend(loc='best',fontsize=size_ticks,frameon=False)
plt.yscale('log')


#Remove frames
sns.despine(top=True, right=True, left=False, bottom=False)

#Save and show plots
Name_of_Figure='L-V_solutions_Snippet'
Extensions_Long=['.pdf','.eps','.png','.svg']
Extensions_Short=['.pdf','.svg']
for ext in Extensions_Short:
    plt.savefig(Output_Path + Name_of_Figure + ext,dpi=300)

plt.show()
#====================================================================
