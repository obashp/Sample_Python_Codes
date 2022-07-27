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
import Utils as Ut

majorLocator = MultipleLocator(0.5)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(0.25)


###----------Amplitude (A10) plot function---------#####

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
 
    
   

folder_path0 = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\D_O'
Ur = GetUrspan(folder_path0)
Urfile = os.path.join(folder_path0,'Urlist.dat')
with open(Urfile, 'r') as f:
    Ur = np.loadtxt(f)
Ur.sort()

Urs = [3.5,4,4.5,5,6,6.5,7,8,10,12,12.5,15,17,19]
Urs1 = [10,12,12.5,15,17,19]
Uri = np.zeros(len(Urs),dtype = int)
nbody = 1

j = 0
for i in range(len(Ur)):
    if Ur[i] == Urs[j]:
        Uri[j] = i;
        j = j+1;
    if j==len(Urs):
        break

rowNum = 2
colNum = 2

figW, figH = myplt.figsize_mm(figSWmm=200, figAR=1.0,
                              nRows=rowNum, nCols=colNum)



fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=rowNum, ncols=colNum,                                                                                                                                                                                                             
                        sharex=False,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=10,
                    wspace=2,hspace=4)
    
x1 = np.zeros((2),dtype = float)
y1 = np.zeros((2),dtype = float)

x1[0] = -3.0
x1[1] = 3.0


k=0
for i in range(1):
    for j in range(1):
        ax1 = axs[1,0]
        ax2 = axs[1,1]
        ax1.plot(y1,x1,'grey')
        ax1.plot(x1,y1,'grey')
        ax2.plot(y1,x1,'grey')
        ax2.plot(x1,y1,'grey')
        
        
        filepath = '\\D_O_' + str(4)
        filename =  'phasep_body_' + str(1) + '.dat'
        phasefile = os.path.join(folder_path0+filepath,filename)
        with open(phasefile,'r') as f:
            data = np.loadtxt(f,skiprows=2)
    
        n = len(data[:,0])
        time = np.zeros(n,dtype=float)
        y = np.zeros(n,dtype=float)
        vy = np.zeros(n,dtype = float)
        ay =np.zeros(n,dtype = float)
        pdata = np.zeros(n,dtype = int)
        
        y = data[150:-150,0]
        vy = data[150:-150,1]
        ay = data[150:-150,2]
        
        k = 0
        i = 1
        while(i < len(ay)):
            if(ay[i-1] < 0  and ay[i] > 0):
                pdata[k] = i
                k = k+1
            i = i+1
        
        pdata = np.trim_zeros(pdata,'b')
        print(len(pdata))
     #   print(pdata)
        
        
        vyp = np.zeros(len(pdata),dtype = float)
        yp = np.zeros(len(pdata),dtype = float)
        
        for i in range(len(pdata)):
            j = pdata[i]
            vyp[i] = vy[j] + (vy[j-1]-vy[j])*(0-ay[j])/(ay[j-1]-ay[j])
            yp[i] = y[j] + (y[j-1]-y[j])*(0-ay[j])/(ay[j-1]-ay[j])
        
        
        ax2.plot(yp,vyp,'ko')
            
                
            
        
        # if nbody ==1:
        #     y = data[:,0]
        #     vy = data[:,1]
            
        # if nbody ==2:
        #     x = data[:,0]
        #     y = data[:,1]
        #     vy = data[:,3]
        
        
        ax1.plot(y,vy,'k')
        
        filename =  'phasep_body_' + str(2) + '.dat'
        phasefile = os.path.join(folder_path0+filepath,filename)
        with open(phasefile,'r') as f:
            data = np.loadtxt(f,skiprows=2)
    
        n = len(data[:,0])
        x = np.zeros(n,dtype=float)
        vx = np.zeros(n,dtype = float)
        ax = np.zeros(n,dtype = float)
        y = np.zeros(n,dtype=float)
        vy = np.zeros(n,dtype = float)
        ay = np.zeros(n,dtype = float)
        pdata = np.zeros(n,dtype = int)
        
        time = data[:,0]
        x = data[150:-150,0]
        x = x - np.mean(x)
        y = data[150:-150,1]
        vx = data[150:-150,2]
        vy = data[150:-150,3]
        ax = data[150:-150,4]
        ay = data[150:-150,5]
        
        ax1.plot(y,vy,'b')
        
        
        k = 0
        i = 1
        while(i < len(ay)):
            if(ay[i-1] < 0 and ay[i]>0):
                pdata[k] = i+1
                k = k+1
            i = i+1
        
        pdata = np.trim_zeros(pdata,'b')
        print(len(pdata))
        print(pdata)
        
        
        vyp = np.zeros(len(pdata),dtype = float)
        yp = np.zeros(len(pdata),dtype = float)
        
        for i in range(len(pdata)):
            j = pdata[i]
            vyp[i] = vy[j] + (vy[j-1]-vy[j])*(0-ay[j])/(ay[j-1]-ay[j])
            yp[i] = y[j] + (y[j-1]-y[j])*(0-ay[j])/(ay[j-1]-ay[j])
        
        
        ax2.plot(yp,vyp,'bo')
        
        
        
        title= 'Ur_' + str(3.5)
        plotparams(ax1,'$v_y$','$y$',[-1.5,1.5],[-1.5,1.5],True,True)
        plotparams(ax2,'$v_y$','$y$',[-1.5,1.5],[-1.5,1.5],True,False)
        
        ax1 = axs[0,0]
        ax2 = axs[0,1]
        ax1.plot(y1,x1,'grey')
        ax1.plot(x1,y1,'grey')
        ax2.plot(y1,x1,'grey')
        ax2.plot(x1,y1,'grey')
        
        ax1.plot(x,vx,'b')
        
        
        
        pdata = np.zeros(n,dtype = int)
        k = 0
        i = 1
        while(i < len(ax)):
            if(ax[i-1] < 0 and ax[i]>0):
                pdata[k] = i+1
                k = k+1
            i = i+1
        
        pdata = np.trim_zeros(pdata,'b')
        print(len(pdata))
        
        
        vyp = np.zeros(len(pdata),dtype = float)
        yp = np.zeros(len(pdata),dtype = float)
        
        for i in range(len(pdata)):
            j = pdata[i]
            vyp[i] = vx[j] + (vx[j-1]-vx[j])*(0-ax[j])/(ax[j-1]-ax[j])
            yp[i] = x[j] + (x[j-1]-x[j])*(0-ax[j])/(ax[j-1]-ax[j])
        
        
        ax2.plot(yp,vyp,'bo')
        plotparams(ax1,'$v_y$','$y$',x1,x1,True,True,title)
        plotparams(ax2,'$v_y$','$y$',x1,x1,True,False,title)
        
        #plt.plot(x,y,'b')
        
        
        # ax = axs[0,1]
        # ax.plot(x,vx,'k')
        # title= 'Ur_' + str(Urs[k])
        # plotparams(ax,'$vx$','$x$',x1,x1,True,True,title)
        # fig.tight_layout()
        # fig.savefig('Phaseplots_D.jpg', dpi=150,format='jpg',bbox_inches='tight', pad_inches=0.1)
fig.tight_layout()
fig.savefig('Phaseplots_D.jpg', dpi=150,format='jpg',bbox_inches='tight', pad_inches=0.1)
    

    