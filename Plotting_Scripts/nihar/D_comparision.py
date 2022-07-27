import os
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import gridspec
import matplotlib.lines as mlines
import Atul_myplt as myplt

majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1)


###----------Amplitude (A10) plot function---------#####

def plotparams(plotloc,ylabel,xlabel,xlim,ylim,isx):
    plotloc.xaxis.set_minor_locator(minorLocator)
    plotloc.xaxis.set_major_locator(majorLocator)

    plotloc.set_ylabel(ylabel,fontsize=12,color='black')
    bot = 0
    if isx:
        plotloc.set_xlabel(xlabel,fontsize=12,color='black')
        bot=1
        
    plotloc.set_xlim(xlim)
    plotloc.set_ylim(ylim)
    
    plotloc.tick_params(which='both',left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=bot, width=0.5, direction="in")
    


figW, figH = myplt.figsize_mm(figSWmm=280, figAR=0.4,
                              nRows=1, nCols=1)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=1, ncols=1,
                        sharex=True,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=10,
                  wspace=2,hspace=4)


    

folder_path0 = 'F:\\Nihar\\VIV_tests\\D_O'
Ampfile = os.path.join(folder_path0,'D_A_Sty_Ur.dat')
Ur0 = np.zeros(38,dtype=float)
AmaxDO = np.zeros(38,dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur0 = data[:,0]
    AmaxDO[:] = data[:,2]


folder_path0 = 'F:\\Nihar\\VIV_tests\\IsolatedD'
Ampfile = os.path.join(folder_path0,'1y_A_Sty_Ur.dat')   
Ur1 = np.zeros(15,dtype=float)
AmaxD = np.zeros(15,dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur1 = data[:,0]
    AmaxD[:] = data[:,1]
    
    
labels = ['(a)','(b)']
for t in range(2):
    ax = axs[0,t]
    ax.text(-0.05, 0.93, labels[t],
        verticalalignment='top', horizontalalignment='right',   
        transform=ax.transAxes,
        color='black', fontsize=15)
    
    
plotformat = ['bv-.','r+--','kv-','c+--' ]

xlim = [Ur0[0]-0.5,Ur0[len(Ur0)-1]+0.5]
ax = axs[0,0]
ax.plot(Ur0,AmaxDO,plotformat[2])
ax.plot(Ur1,AmaxD,plotformat[1])
#ax.legend(['$DO - O-A_{10}$ (Present study)','$OO - O-A_{10}$ (Prasanth et.al)','$OO-O-A_{10}$ (Bao et.al)','$O-A_{10}$ (Bao et.al)'])
plotparams(ax,'$A_{10}$','U_R',xlim,[0.0,4.0],True)
# ax.plot(Ur,phasepi,linestyle='--',color='grey')  
# ax.plot(Ur,phasepi*0.5,linestyle='--',color='grey')  
# ax.plot(Ur,phasepi*0,linestyle='--',color='grey')  



     
fig.tight_layout()
fig.savefig('AmpcompD.jpg', dpi=100,format='jpg',bbox_inches='tight', pad_inches=0.1)