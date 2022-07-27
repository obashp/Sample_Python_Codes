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


def plotparams(plotloc,ylabel,xlabel,xlim,ylim,isx=True,isy=True,plttitle=''):
    plotloc.xaxis.set_minor_locator(minorLocator)
    plotloc.xaxis.set_major_locator(majorLocator)

    plotloc.title.set_text(plttitle)
    bot = 0
    let = 0
    rit = 0
    if isx:
        plotloc.set_xlabel(xlabel,fontsize=12,color='black')
        bot=1
    if isy:
        plotloc.set_ylabel(ylabel,fontsize=12,color='black')
        let=1
    else:
        rit=1
        
    plotloc.set_xlim(xlim)
    plotloc.set_ylim(ylim)
    
    plotloc.tick_params(which='both',left=1, top=1, right=1,
                    bottom=1, labelleft=let, labeltop=0,
                    labelright=rit, labelbottom=bot, width=0.5, direction="in")


figW, figH = myplt.figsize_mm(figSWmm=260, figAR=0.4,
                              nRows=2, nCols=1)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=2, ncols=1,
                        sharex=True,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=10,
                  wspace=2,hspace=4)

folder_path0 = 'F:\\Nihar\\VIV_tests\\grid_resolution'
plotformat = ['k-.','r-.','g:','c--' ]
gridsize= ['$257x193 (0.04)$','$385x257 (0.02)$','$449x257 (0.0175)$','$513x321 (0.015)$']
lwt=3
xlim = [440,470]
ylim = [-2,5]
ylim1 = [-2,4]
ins = 0
ins1 = 0

for t in range(4):
    fnameD = 'data_' + str(t+1) + '.dat'
    fnameO = 'data_' + str(t+5) + '.dat'
    
    AmpfileD = os.path.join(folder_path0,fnameD)
    with open(AmpfileD,'r') as f:
        data = np.loadtxt(f,skiprows=1)
        
    
    ax = axs[0,0]
    ax.plot(data[:,0],data[:,1],plotformat[t],linewidth=lwt)
   
    
    plotparams(ax,'y','time',xlim,ylim,False,True,'D-Cylinder')
    
    if ins==0:
        ins = ax.inset_axes([0.3,0.65,0.25,0.3])
        ins.set_xlim([452,453.7])
        ins.set_ylim([1.25,1.65])
        ax.indicate_inset_zoom(ins)
    
    ins.plot(data[:,0],data[:,1],plotformat[t],linewidth=lwt)
    
    
        
    AmpfileO = os.path.join(folder_path0,fnameO)
    with open(AmpfileO,'r') as f:
        data = np.loadtxt(f,skiprows=1)
        ax = axs[1,0]
        ax.plot(data[:,0],data[:,1],plotformat[t],linewidth=lwt)
        
    
    plotparams(ax,'y','time',xlim,ylim1,True,True,'O-Cylinder')
    
    if ins1==0:
        ins1 = ax.inset_axes([0.3,0.65,0.25,0.3])
        ins1.set_xlim([460,461.7])
        ins1.set_ylim([0.75,1.1])
        ax.indicate_inset_zoom(ins1)
        ax.title.set_text('O-Cylinder')
    
    ins1.plot(data[:,0],data[:,1],plotformat[t],linewidth=lwt)
    
fig.legend(gridsize,framealpha=1.0)
fig.tight_layout()
fig.savefig('Grid_Resolution.jpg', dpi=150,format='jpg',bbox_inches='tight', pad_inches=0.1)