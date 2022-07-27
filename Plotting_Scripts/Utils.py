import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

def safe_ln(x, minval=0.001):
    return np.log(x.clip(min=minval))

majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1.0)

def plotparams(plotloc,ylabel,xlabel,xlim,ylim,isx = False, yformatter = False, yformval = 0.0, isy= True, xformatter = False, xformval = 0.0):
    
    if(xformatter):
        xmajorLocator = MultipleLocator(xformval)
        xmajorFormatter = FormatStrFormatter('%d')
        xminorLocator = MultipleLocator(0.5*xformval)
        plotloc.xaxis.set_minor_locator(xminorLocator)
        plotloc.xaxis.set_major_locator(xmajorLocator)
    else:
        plotloc.xaxis.set_minor_locator(minorLocator)
        plotloc.xaxis.set_major_locator(majorLocator)
    
    if(yformatter):
        ymajorLocator = MultipleLocator(yformval)
        ymajorFormatter = FormatStrFormatter('%d')
        yminorLocator = MultipleLocator(0.2*yformval)
        plotloc.yaxis.set_minor_locator(yminorLocator)
        plotloc.yaxis.set_major_locator(ymajorLocator)

    
    if isy:
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
    
def GetUrspan(folder_path):
    Urfile = os.path.join(folder_path,'Urlist.dat')
    with open(Urfile, 'r') as f:
       Ur = np.loadtxt(f)   
       
    Ur.sort()
    return Ur

def plotaxes(plotloc, x1, x2, y1, y2):
    
    plotloc.plot([0.5*(x1+x2),0.5*(x1+x2)],[y1,y2],'grey')
    plotloc.plot([x1,x2],[0.5*(y1+y2),0.5*(y1+y2)],'grey')
    
def TimeSeriesplot(plotloc,file,tidx,yidx,lstyle):
    with open(file,'r') as f:
        data = np.loadtxt(f,skiprows=1)
        time = data[:,tidx]
        y = data[:,yidx]
        plotloc.plot(time,y,lstyle)
    

def Amplitudeplot(plotloc, file, nx, ft, plotlabel,isx = False,yformatter = False,yformval = 0.0):
    
    AmaxA = np.zeros(nx,dtype=float)
    AmaxM = np.zeros(nx,dtype=float)
    Amaxm = np.zeros(nx,dtype=float)
    error = np.zeros((2,nx),dtype=float)
    with open(file,'r') as f:
        data = np.loadtxt(f,skiprows = 1)
        Ur = data[0:nx,0]
        AmaxA[:] = data[0:nx,2]
        AmaxM[:] = data[0:nx,1]
        Amaxm[:] = data[0:nx,3]
        #print(Amaxm[-1],Ur[-1],nx)
        xlim = [Ur[0],Ur[-1]]
        ylim = [0.0,max(AmaxA)+0.3]
        error[1,:] = AmaxM[:] - AmaxA[:]
        error[0,:] = -Amaxm[:] + AmaxA[:]
        plotloc.errorbar(Ur,AmaxA,yerr=error,fmt=ft,ecolor='k',capsize=3.5)
        plotparams(plotloc,'$A_{10}$','U_R',xlim,ylim,isx,yformatter,yformval)
        plotloc.text(-0.05, 0.93, plotlabel,
        verticalalignment='top', horizontalalignment='right',
        transform=plotloc.transAxes,
        color='black', fontsize=12)
    ratio = (ylim[1]-ylim[0])/(xlim[1]-xlim[0])
    return ratio
        
def Cfplot(plotloc, file, nx, ft, nvar, phaseflag = False, isx = False):
    
    y = np.zeros(nx,dtype=float)
    with open(file, 'r') as f:
        data = np.loadtxt(f,skiprows = 1)
        Ur = data[:,0]
        y = data[:,nvar]
        if(phaseflag):
            for i in range(len(y)):
                if(y[i] > 180):
                    y[i] = 360-y[i]
                elif(y[i] < 0):
                    y[i] = -y[i]
        plotloc.plot(Ur,y,ft, markersize = 4.2)
    return max(y)
    
def plotconst(plotloc, x, const):
    cvec = np.zeros(len(x),dtype = float)
    for i in range(len(x)):
        cvec[i] = const
    plotloc.plot(x,cvec,linestyle = '--',color = 'grey')
    
def phasedata(filepre, folder_path, folderpre, Ur, xvar, yvar):
    
    filepath = folderpre + str(Ur) + '\\'
    filename = filepre + '_phasep_body.dat'
    file = folder_path + filepath + filename
    with open(file,'r') as f:
        data = np.loadtxt(f,skiprows = 2)
    x = data[3000:-3000,xvar]    
    y = data[3000:-3000,yvar]
    
    return x,y

def equalize(x1,x2):
    
    i1 = 0
    i2 = 0
    #print(x1[0],x2[0])
    if(x1[0] > x2[0]):
        i1 = 0
        for i in range(len(x2)):
            if(x2[i] == x1[0]):
                i2 = i
                break
    
    if(x2[0] >= x1[0]):
        i2 = 0
        for i in range(len(x1)):
            if(x1[i] == x2[0]):
                i1 = i
                break
            
    return i1,i2
            
            

        
def PSDplot(plotloc, file, Ur, St, NUMVAR, ylabel,plotlabel, isx = False, logswitch = True):
    
    NI = len(Ur)
    NJ = 0    
    
    
    
    with open(file,'r') as f:
        while NJ == 0:
            line = f.readline()
            line = line.split()
            if(line[0] == 'zone'):
                NJ = int(line[2])
    #print(NI,NJ)
    xc = np.zeros((NI,1),dtype=float)
    yc = np.zeros((NJ,1),dtype=float)
    qq = np.zeros((NJ,NI),dtype=float)
        
      
    with open(file, 'r') as f: 
        data1 = np.loadtxt(f,skiprows=2)
 
    x = data1[:,0]
    y = data1[:,1]

    yll = y.min();  yul = y.max()
    xc = Ur 
    yc = np.linspace(yll, yul, NJ) 
    X,Y = np.meshgrid(xc,yc)

    levels = [-4,-3,-2.75,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0]
    
    
    count = 0
    for i in range(NI):
        for j in range(NJ):
            qq[j,i] = data1[count,NUMVAR]
            count = count+1
   
    qq[:,:] = np.minimum(qq[:,:],3)
    qq[:,:] = np.maximum(qq[:,:],0)
    fq = qq[:,:]
    if logswitch == True:
        fq = safe_ln(fq)
    
    
    ax1 = plotloc.contourf(X,Y,fq,100,cmap=plt.cm.afmhot_r) 
    #ax1 = plotloc.contourf(X,Y,fq,levels,cmap=plt.cm.afmhot_r)
    
    xlim = [Ur[0],Ur[-1]]
    ylim = [min(y),max(y)]
    
    plotparams(plotloc,ylabel,'$U_R$',xlim,ylim,isx,True,1.0,True,True,1.0)
    plotloc.text(-0.05, 0.93, plotlabel,
        verticalalignment='top', horizontalalignment='right',
        transform=plotloc.transAxes,
        color='black', fontsize=12)
    
    return ax1
    
    
def plotStfN(plotloc, StfNflag, St, Stf, Ur, scalefact = 1.0):
    
    if(Stf):
        plotloc.text(0.64, 0.62, '$f_{v}$ =  ' + str(round(St,3)),
        verticalalignment='top', horizontalalignment='right',
        transform=plotloc.transAxes,
        color='black', fontsize=12)
        
    nUr = len(Ur)
    if(StfNflag):
        St1 = np.zeros(nUr)
        fN = np.zeros(nUr)
        for i in range(nUr):
            St1[i] = scalefact*St*Ur[i]
            fN[i] = 1.0
        plotloc.plot(Ur,St1, linestyle='--', label = '', color = 'k')
        
      
    # if(StfNflag):
    #     plotloc.plot(Ur,fN, linestyle='--', label = '', color = '0.7')
    #     plotloc.plot(Ur,3*fN, linestyle='--', label = '', color = '0.7')
    #     plotloc.plot(Ur,5*fN, linestyle='--', label = '', color = '0.7')
        
def plot2Dgrid(plotloc,xfile,yfile,stylestring):
    stylestring = stylestring.split()
    colour = stylestring[0]
    style = stylestring[1]
    
    with open(xfile,'r') as f:
        xdata = np.loadtxt(f)
        x = xdata[:,1]
    
    with open(yfile,'r') as f:
        ydata = np.loadtxt(f)
        y = ydata[:,1]
    
    nx = len(xdata[:,0])
    ny = len(ydata[:,0])
    print(nx,ny)
    del xdata
    del ydata

    ymax = max(y)
    ymin = min(y)
    print(ymin,ymax)

    xmax = max(x)
    xmin = min(x)    
    print(xmin,xmax)

    for i in range(nx):
        plotloc.plot([x[i],x[i]],[ymin,ymax],color='xkcd:black',linestyle='solid',linewidth='0.2')
    
    for i in range(ny):
        plotloc.plot([xmin,xmax],[y[i],y[i]],color='xkcd:black',linestyle='solid',linewidth='0.2')
        
    plotparams(plotloc,'Y','X',[xmin,xmax],[ymin,ymax],True,True,2.0,True,True,2.0)
    
def plot2Dgridzoomed(plotloc,xfile,yfile,stylestring,stix,edix,stiy,ediy):
    stylestring = stylestring.split()
    colour = stylestring[0]
    style = stylestring[1]
    
    with open(xfile,'r') as f:
        xdata = np.loadtxt(f)
        x = xdata[:,1]
    
    with open(yfile,'r') as f:
        ydata = np.loadtxt(f)
        y = ydata[:,1]
    
    del xdata
    del ydata
    
    xmin = x[stix-1]
    xmax = x[edix-1]
    ymin = y[stiy-1]
    ymax = y[ediy-1]
    
    i = stix-1
    while i < edix:
        plotloc.plot([x[i],x[i]],[ymin,ymax],color='xkcd:black',linestyle='solid',linewidth='0.2')
        i = i+1
        
    i = stiy-1
    while i < ediy:
        plotloc.plot([xmin,xmax],[y[i],y[i]],color='xkcd:black',linestyle='solid',linewidth='0.2')
        i = i+1
        
    
        
    
    
   
def plot2DUnstrucSurf(plotloc,filename,stylestring):
    stylestring = stylestring.split()
    colour = stylestring[0]
    style = stylestring[1]
    
    N = 0
    with open(filename,'r') as f:
        while N == 0:
            line = f.readline()
            line = line.split()
            if(line[0] == 'zone'):
                N = int(line[2])
        
        line = f.readline()
        N = int(N/3)
        print(N)
        x = np.zeros(N,float)
        y = np.zeros(N,float)
        i = 0
        while i < N:
            line = f.readline()
            line = line.split()
            x[i] = line[0]
            y[i] = line[1]
            i = i+1

        plotloc.plot(x,y,color = colour,linestyle=style,linewidth = 0.5)
        

    
    
    

