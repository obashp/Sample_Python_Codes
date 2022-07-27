# -*- coding: utf-8 -*-
import os
import numpy as np
import numpy.ma as ma
import pylab as pl
import matplotlib.pyplot as plt

import matplotlib.patches as mpatches
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import gridspec
import matplotlib.lines as mlines
import Atul_myplt as myplt
import Atul_read_vicar_data as rvd
import streams as st

folder_path0 = 'I:\\NIHAR\\Tandem_D_O\\D_O\\D_O_6.5\\plt'
folder_path1 = 'I:\\NIHAR\\Tandem_D_O\\D_O\\D_O_6.5'

gridinfo = os.path.join(folder_path1,'gscript.dat')
with open(gridinfo,'r') as f:
    data = np.genfromtxt(f,delimiter=',',skip_header=1,skip_footer=2)
    
    
pltfilelist = ['q0075200.plt','q0075300.plt','q0075400.plt','q0075500.plt']
markerfilelist = ['marker.0075200.dat','marker.0075300.dat','marker.0075400.dat','marker.0075500.dat'] 

figW, figH = myplt.figsize_mm(figSWmm=125, figAR=1/2,
                              nRows=4, nCols=2)



fig, axs = plt.subplots(figsize=(figW, figH),nrows=4, ncols=2, 
                        gridspec_kw={'width_ratios': [1, 3],'height_ratios': [1, 1, 1, 1]},sharey=True)    
#figsize=(figW, figH),
#,gridspec_kw ={'width_ratios': [4, 7]}
myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=8,
                    wspace=4,hspace=4)

for i in range(len(pltfilelist)):
    
    
    
    A1 = os.path.join(folder_path0,pltfilelist[i])
    A2 = os.path.join(folder_path0,markerfilelist[i])
    
    xc,yc,qqf = rvd.read_binary(A1)
    
    #plt variables = xc,yc,u,v,p,div,vor,bl,gc
    
     
    X,Y = np.meshgrid(xc,yc)
    
    qqf[:,:,0,4] = np.minimum(qqf[:,:,0,4],2.0)
    qqf[:,:,0,4] = np.maximum(qqf[:,:,0,4],-2.0)
    
    mask=(qqf[:,:,0,5] > 0.5)
    vort = ma.masked_array(qqf[:,:,0,4],mask = mask)
            
    ax = axs[i,1]
    ax.contourf(X,Y,vort,50,cmap=plt.cm.bwr)
    xlim = [4.5,20.5]
    ylim = [6.5,13.5]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ar = (ylim[1]-ylim[0])/(xlim[1]-xlim[0])
    ax.set_aspect(1)
    ax.set_adjustable("box")
    ax = axs[i,0]
    
    qqf[:,:,0,2] = np.minimum(qqf[:,:,0,2],0.5)
    qqf[:,:,0,2] = np.maximum(qqf[:,:,0,2],-0.5)
    
    pres = ma.masked_array(qqf[:,:,0,2],mask = mask)
    ax.contourf(X,Y,pres,50,cmap=plt.cm.bwr)
    xlim = [4.5,12.5]
    ylim = [6.5,13.5]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ar = (ylim[1]-ylim[0])/(xlim[1]-xlim[0])
    ax.set_aspect(1)
    ax.set_adjustable("box")
    
    # rect = [i+1, i, i]    
    # divider = HBoxDivider(
    #     fig, rect,
    #     horizontal=[Size.AxesX(axs[i,0]), Size.Fixed(0.5), Size.AxesX(axs[i,1])],
    #     vertical=[Size.AxesY(axs[i,0]), Size.Scaled(1), Size.AxesY(axs[i,1])])
    # ax[i,0].set_axes_locator(divider.new_locator(0))
    # ax[i,0].set_axes_locator(divider.new_locator(2))
    
    
    lx = len(xc)
    ly = len(yc)
    lx1 = lx//6
    ly1 = ly//6
    xp = np.zeros(lx1)
    xind = np.zeros(lx1,dtype = int)
    yp = np.zeros(ly1)
    yind = np.zeros(ly1,dtype = int)
    
    
    
    i = 0 
    j=0
    while i < len(xc) and j < lx1:
        xp[j] = xc[i]
        xind[j] = i
        j = j+1
        i = i+6
        
    i = 0
    j = 0
    while i < len(yc) and j < ly1:
        yp[j] = yc[i]
        yind[j] = i
        j = j+1
        i = i+6
        
    u = np.zeros(lx1*ly1)
    v = np.zeros(lx1*ly1)
    
    k =0
    for t in xind:
        for p in yind:
            u[k] = qqf[p,t,0,0]
            v[k] = qqf[p,t,0,1]
            k = k+1
            
    Xp,Yp = np.meshgrid(xp,yp)
    ax.quiver(Xp,Yp,u,v)
    
if os.path.exists("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Vort_6.5.pdf"):
  os.remove("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Vort_6.5.pdf")

fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Vort_6.5.pdf', dpi=300,format='pdf',bbox_inches='tight', pad_inches=0.1)



