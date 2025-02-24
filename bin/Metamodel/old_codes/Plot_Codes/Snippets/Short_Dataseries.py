#2022-01-17. A snippet to plot short datasets. For instance, you simulate a system for 4/values of a parameter


import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import figure 

#Parameters and fontsizes
#========================
size_axis=14
size_ticks=12
size_title=16 
#========================

#Manually introduced data
#================================================================================
Series1=[0.723,0.741,0.768,0.772]#S-8
Series2=[0.719,0.739,0.760,0.763]#V-10
Series3=[0.717,0.722,0.728,0.729]#V-22
Series4=[0.649,0.660,0.681,0.685]#V-23_24
Series5=[0.720,0.731,0.750,0.753]#V-25

Threshold=[1e-4,5e-5,1e-5,5e-6]
#================================================================================

#Plot Figure
#================================================================================
figure(figsize=(9, 6), dpi=80)

plt.plot(Threshold,Series1,'-o',color='blue',linewidth=3,label='Label1')
plt.plot(Threshold,Series2,'-o',color='orange',linewidth=3,label='Label2')
plt.plot(Threshold,Series3,'-o',color='green',linewidth=3,label='Label3')
plt.plot(Threshold,Series4,'-o',color='red',linewidth=3,label='Label4')
plt.plot(Threshold,Series5,'-o',color='purple',linewidth=3,label='Label5')

#Axes, title and ticks
plt.ylabel('y-axis-name',fontsize=size_axis)
plt.xlabel('x-axis-name',fontsize=size_axis)


plt.xticks(fontsize=size_ticks);plt.yticks(fontsize=size_ticks)
xticks_pos=Threshold
xticks_labels=Threshold

ymin=0.62;ymax=0.80
#Ylimitis
plt.ylim(ymax = ymax, ymin = ymin)

#Legend and scale
plt.xscale('log')
plt.legend(loc='upper right')

#Save plots in png and pdf
#plt.savefig(Output_Path+'Name.png',dpi=300)
#plt.savefig('Name.pdf',dpi=300)

#axes limits
#plt.ylim(-2, 2)
#plt.xlim(0,10)

plt.show()
#==============================================================================
