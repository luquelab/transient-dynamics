import numpy as np
import matplotlib
from decimal import Decimal
from matplotlib.pyplot import figure
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import sys

Min_val = sys.argv[1]
Max_val = sys.argv[2]
cmap=sys.argv[3]

Min_val=float(Min_val)
Max_val=float(Max_val)
print Min_val
print Max_val
path_plot='/home/sergio/work/Github/needle-finder/results/Plots/Paper_Matt/'

#Min and Max Values for colorbar
# Dataframe=pd.read_csv(path_plot+'Data_colorbar.csv', index_col=0)
# print Dataframe
# print Dataframe.drop
# # Dataframe.to_csv=(path+'Data_colorbar.csv')  
# Dataframe_transp = Dataframe.transpose()
# print Dataframe_transp
# #Viral subdataset                               
# Dataframe_virus=Dataframe_transp.iloc[:2]
# Max_val_0_virus=Dataframe_virus.max()               
# Max_val_virus=np.max(Dataframe_virus.max())    
# Min_val_0_virus=Dataframe_virus.min()                             
# Min_val_virus=np.min(Dataframe_virus.min())
# #Bacterial subdataset
# Dataframe_bacteria=Dataframe_transp.iloc[2:] 
# Max_val_0_bacteria=Dataframe_bacteria.max()  
# Max_val_bacteria=np.max(Dataframe_bacteria.max())   
# Min_val_0_bacteria=Dataframe_bacteria.min()               
# Min_val_bacteria=np.min(Dataframe_bacteria.min())

# Max_val=max(Max_val_bacteria, Max_val_virus)
# Min_val=max(Min_val_bacteria, Min_val_virus)
# # Max_val_0=Dataframe.max()
# # Max_val=np.max(Dataframe.max())
# # Min_val_0=Dataframe.min()
# # Min_val=np.min(Dataframe.min())


size_ticks=5

cm = 1/2.54
figure(figsize=(0.4*cm, 1.6*cm), dpi=300)
colorbar=plt.subplot()
colbar_ticks=[0,0.5,1]


cb1=matplotlib.colorbar.ColorbarBase(colorbar, cmap=cmap,ticks=(colbar_ticks))
cb1.set_ticks(colbar_ticks)

#Ticks for minimum value being very close to zero

#cb1.set_ticklabels(("{:.2f}".format(Decimal(Min_val)),"{:.2f}".format(0.5*Max_val),"{:.2f}".format(Max_val)))

cb1.set_ticklabels(("{:.0f}".format(Decimal(Min_val)),"{:.0f}".format(0.5*Max_val),"{:.0f}".format(Max_val)))

#cb1.set_ticklabels((0,1,2))
cb1.ax.tick_params(labelsize=size_ticks)
cb1.outline.set_visible(False)


Extensions=['.pdf','.eps','.png','.svg']
Name_of_Figure='Colorbars_' + str(cmap)

for ext in Extensions:
    plt.savefig(path_plot+Name_of_Figure + ext,dpi=300,bbox_inches='tight',transparent=True)
plt.show()

cb1.ax.tick_params(labelsize=size_ticks)
cb1.outline.set_visible(False)
