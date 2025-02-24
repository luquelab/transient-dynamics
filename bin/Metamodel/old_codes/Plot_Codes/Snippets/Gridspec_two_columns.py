#2021-02-02. This code only works with python 2.


#Import libraries
#+++++++++++++++++++++++++++++++++++++++++++++++++++
#from Metamodel_Functions import *
import numpy as np
from random import shuffle
import matplotlib
import matplotlib.gridspec as gridspec
import time
import matplotlib.pyplot as plt
from decimal import Decimal
import seaborn as sns
#+++++++++++++++++++++++++++++++++++++++++++++++++++


#PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0. Fontsizes
#========================
size_axis=14
size_ticks=11
size_letter=15
#=======================

#1. Configure gridspec plot
#===============================================
cm=1/2.54# centimeters in inches
Rows=1;Cols=2
fig=plt.figure(figsize=(8*cm, 4*cm))
gs=gridspec.GridSpec(Rows,Cols, width_ratios=[1,0.04])
gs.update(left=0.07,right=0.85,bottom=0.08,top=0.97,wspace=0.12,hspace=0.02)
#===============================================

#2. Heatmap (0,0)
#===============================================
ax_00=plt.subplot(gs[0,0])
uniform_data = np.random.rand(10, 12)
sns.heatmap(uniform_data,vmin=0,vmax=1,cmap="Blues",cbar=False)
#===============================================  

#3.4. (1,0) - Colorbar
#========================================================================
colorbar=plt.subplot(gs[0,1])

colbar_ticks=[0,0.5,1]
colmap='Blues'
cb1=matplotlib.colorbar.ColorbarBase(colorbar, cmap=colmap,ticks=(colbar_ticks)) 
cb1.set_ticks(colbar_ticks)

cb1.set_ticklabels(("{:.2f}".format(Decimal(0)),"{:.2f}".format(0.5),"{:.2f}".format(1))) 

cb1.ax.tick_params(labelsize=size_ticks)
cb1.outline.set_visible(False)
cb1.set_label("Title",fontsize=size_axis)
#===============================================

plt.show()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


