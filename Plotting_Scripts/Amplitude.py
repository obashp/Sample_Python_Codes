# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:23:03 2020

@author: Laxman
"""


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
    



figW, figH = myplt.figsize_mm(figSWmm=260, figAR=0.4,
                              nRows=3, nCols=1)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=3, ncols=1,
                        sharex=True,squeeze=False)  # sharing the x-label

fig1, axs1 = plt.subplots(figsize=(figW, figH),
                        nrows=3, ncols=1,
                        sharex=True,squeeze=False)

myplt.margins_adjust(fig, figW, figH, left=10, right=10, top=5, bottom=5,
                   wspace=0,hspace=2)
myplt.margins_adjust(fig1, figW, figH, left=10, right=10, top=5, bottom=5,
                   wspace=0,hspace=2)

folder_path = ['E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\D_O']

plotformat = ['bv-.','r+--','kv-','c+--' ]
print(folder_path)

pindex = 0
while pindex < len(folder_path):
    Ampfilelist = os.path.join(folder_path[pindex],'Afilelist.dat')
    with open(Ampfilelist,'r') as f:
        filelist = np.loadtxt(f,dtype=str) 
    nBodies = int(filelist[0])
    dof = np.zeros(nBodies,dtype=int)
    iBody = 1
    while iBody <= nBodies:
        dof[iBody-1] = filelist[iBody]
        iBody = iBody+1

    nUr = int(filelist[-1])
    filelist = filelist[iBody:-1]
    nfiles = len(filelist) 

    AmaxA = np.zeros((nfiles,nUr),dtype=float)
    AmaxM = np.zeros((nfiles,nUr),dtype=float)
    Amaxm = np.zeros((nfiles,nUr),dtype=float)
    Cxrms = np.zeros((nfiles,nUr),dtype=float)
    Cyrms = np.zeros((nfiles,nUr),dtype=float)
    Cvrms = np.zeros((nfiles,nUr),dtype=float)
    phasey = np.zeros((nfiles,nUr),dtype=float)
    phasex = np.zeros((nfiles,nUr),dtype=float)
    phasev = np.zeros((nfiles,nUr),dtype=float)
    phasepi = np.zeros(nUr,dtype=float)
    for i in range(nUr):
        phasepi[i] = 180

    yc = 0
    xc = 0
    for iName in range(nfiles):
        Ampfile = os.path.join(folder_path[pindex],filelist[iName])
        body = int(filelist[iName][0])
        dcode = dof[body-1]
        dirn = filelist[iName][1]
        # print(dcode,dirn)
        Cyflag = False
        Cxflag = False
        with open(Ampfile,'r') as f:
            data = np.loadtxt(f,skiprows=1)
            Ur = data[:,0]
            print(len(Ur))
            xlim = [Ur[0]-0.5,Ur[len(Ur)-1]+0.5]
            AmaxA[iName,:] = data[:,1]
            AmaxM[iName,:] = data[:,2]
            Amaxm[iName,:] = data[:,3]
            
            if((dcode == 2 or dcode == 12) and dirn == 'y'):
                phasey[iName,:] = data[:,5]
                Cyrms[iName,:] = data[:,6]
                phasev[iName,:] = data[:,7]
                Cvrms[iName,:] = data[:,8]
                Cyflag = True
                if(dcode == 2):
                    Cxrms[iName,:] = data[:,9]
                    Cxflag = True
                
            if(dcode == 12 and dirn == 'x'):
                phasex[iName,:] = data[:,5]
                Cxrms[iName,:] = data[:,6]
                Cxflag = True
                
                
            error = np.zeros((2,nUr),dtype=float)
            error[1,:] = AmaxM[iName,:] - AmaxA[iName,:]
            error[0,:] = -Amaxm[iName,:] + AmaxA[iName,:]
            
            ax = axs[0,0]
            ax1 = axs1[0,0]
            if(dirn == 'y'):
                ax.errorbar(Ur,AmaxA[iName],yerr=error,fmt=plotformat[yc],ecolor='k',capsize=3.5)
                ax.legend(['$D-A_{y}$','$O-A_{y}$'])
                plotparams(ax,'$A_{10}$','U_R',xlim,[0.0,1.5],False)
            elif(dirn == 'x'):
                ax1.errorbar(Ur,AmaxA[iName],yerr=error,fmt=plotformat[xc],ecolor='k',capsize=3.5)
                ax1.legend(['$O-A_{x}$'])
                plotparams(ax1,'$A_{10}$','U_R',xlim,[0.0,0.5],False)
               
             
            ax = axs[1,0]
            if(Cyflag):
                ax.plot(Ur,Cyrms[iName],plotformat[yc])
                ax.plot(Ur,Cvrms[iName],plotformat[yc+2])
                ax.legend(['$D-Cy_{RMS}$','$D-Cv_{RMS}$','$O-Cy_{RMS}$','$O-Cv_{RMS}$'])
                plotparams(ax,'Cy','U_R',xlim,[0.0,1.75],False)
            
            
            if(Cxflag):
                ax1 = axs1[1,0]
                ax1.plot(Ur,Cxrms[iName],plotformat[xc])        
                ax1.legend(['$D-Cx_{RMS}$','$O-Cx_{RMS}$'])
                plotparams(ax1,'Cx','U_R',xlim,[0.0,5],False)
                xc = xc+1
            
            ax = axs[2,0]
            if(dirn == 'y'):
                ax.plot(Ur,phasey[iName],plotformat[yc])
                ax.plot(Ur,phasev[iName],plotformat[yc+2])
                ax.legend(['$D-\u03A6_{yC}$','$D-\u03A6_{yCv}$','$O-\u03A6_{yC}$','$O-\u03A6_{yCv}$'])
                plotparams(ax,'$\u03A6_{yC}$','U_R',xlim,[-30,210],True)
                yc = yc+1
            
            ax1= axs1[2,0]
            if(dirn == 'x'):
                ax1.plot(Ur,phasex[iName],plotformat[xc-1])
                ax1.legend(['$O-\u03A6_{xC}$'])
                plotparams(ax1,'$\u03A6_{xC}$','U_R',xlim,[-30,210],True)
    pindex = pindex+1
       
        
  
ax.plot(Ur,phasepi,linestyle='--',color='grey')  
ax.plot(Ur,phasepi*0.5,linestyle='--',color='grey')  
ax.plot(Ur,phasepi*0,linestyle='--',color='grey')


labels = ['(a)','(b)','(c)']
for t in range(3):
    ax = axs[t,0]
    ax.text(-0.05, 0.93, labels[t],
        verticalalignment='top', horizontalalignment='right',   
        transform=ax.transAxes,
        color='black', fontsize=15)
    ax = axs1[t,0]
    ax.text(-0.05, 0.93, labels[t],
        verticalalignment='top', horizontalalignment='right',   
        transform=ax.transAxes,
        color='black', fontsize=15)

ax1.plot(Ur,phasepi,linestyle='--',color='grey')  
ax1.plot(Ur,phasepi*0.5,linestyle='--',color='grey')  
ax1.plot(Ur,phasepi*0,linestyle='--',color='grey')  



     
# fig.tight_layout()
fig.savefig('Amp_y.pdf', dpi=100,format='pdf',bbox_inches='tight', pad_inches=0.1)
fig1.savefig('Amp_x.pdf', dpi=100,format='pdf',bbox_inches='tight', pad_inches=0.1)



# # Ampfile = os.path.join(folder_path0,'Phdiff_DO_Cy_Cv.dat')
# # with open(Ampfile,'r') as f:
# #     data = np.loadtxt(f,skiprows=1)
# #     phase = data[:,1]
# #     phasev = data[:,2]
    
# # ax.plot(Ur,phase,'bx')
# # ax.plot(Ur,phasev,'ko')
# # ax.legend(['$D-\u03A6_{CyDO}$','$D-\u03A6_{CvDO}$'])
# # plotparams(ax,'$\u03A6_{CDO}$','U_R',xlim,[-90,270],True)  
# # ax.plot(Ur,phasepi,linestyle='--',color='grey')  
# # ax.plot(Ur,phasepi*0.5,linestyle='--',color='grey')  
# # ax.plot(Ur,phasepi*0,linestyle='--',color='grey')   

    
    