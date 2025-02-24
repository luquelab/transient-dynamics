#2022-01-28. A snippet to plot functions. For instance, plot Unique  k-mer metagenome vs read length and kmer size

import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import figure
import numpy as np

#Manually introduced data
#================================================================================
Read_Length=[50,100,200]
kmersize=[np.arange(2,R_L) for R_L in Read_Length]

Hypothetical_Unique_Metagenome=[2*k*(R_L - k + 1) for (R_L,k) in zip(Read_Length,kmersize)]
#================================================================================

#Plot Figure
#====================================================================

#Figure size and resolution, text sizes
#......................................................
cm = 1/2.54  # centimeters in inches
Width=8*cm;Height=4*cm #Width and height of plots
figure(figsize=(Width, Height), dpi=300) #Resolution

size_axis=7;size_ticks=5;size_title=5
#......................................................

#Axes, title and ticks
#......................................................
plt.ylabel('Hypothetical_Unique_Metagenome',fontsize=size_axis)
plt.xlabel('K-mer size',fontsize=size_axis)
#......................................................

#......................................................
plt.plot(kmersize[0],Hypothetical_Unique_Metagenome[0],color='red',label='1',linewidth=3)
plt.plot(kmersize[1],Hypothetical_Unique_Metagenome[1],color='blue',label='2',linewidth=3)
plt.plot(kmersize[2],Hypothetical_Unique_Metagenome[2],color='green',label='3',linewidth=3)
#......................................................

#Log scale
#plt.yscale('log')

#axes limits
y_minimum=0;y_maximum=1.1*np.max(Hypothetical_Unique_Metagenome[2])
plt.ylim(y_minimum, y_maximum)

#Despine figure and legend
plt.legend(loc='upper right', frameon=False)
sns.despine(top=True, right=True, left=False, bottom=False)

#Save figure in all extensions
Name_of_Figure='Name'
Extensions=['.pdf','.eps','.png','.svg']
for ext in Extensions:
    plt.savefig(Output_Path+Name_of_Figure + ext,dpi=300)

plt.show()
#================================================================================
