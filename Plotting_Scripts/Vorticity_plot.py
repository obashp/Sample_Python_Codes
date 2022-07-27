# -*- coding: utf-8 -*-
import os
import numpy as np
import numpy.ma as ma
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib

import matplotlib.patches as mpatches
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import gridspec
import matplotlib.lines as mlines
import Atul_myplt as myplt
import Atul_read_vicar_data as rvd
import Utils as ut
import streams as st


folder_path = 'I:\\NIHAR\\Tandem_D_O\\D_O\\D_O_3.5'
dispdatafile = ['viv_body1dir2_probe.dat','viv_body2dir2_probe.dat']
liftdatafile = ['drag_lift_body_001.dat','drag_lift_body_002.dat']

DF0 = os.path.join(folder_path,dispdatafile[0])
LF0 = os.path.join(folder_path,liftdatafile[0])
DF1 = os.path.join(folder_path,dispdatafile[1])
LF1 = os.path.join(folder_path,liftdatafile[1])

with open(DF0,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    t = data[:,1]
    yD = data[:,2]
    
with open(LF0,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    cyD = data[:,6]
    
with open(DF1,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    yO = data[:,2]

with open(LF1,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    cyO = data[:,6]
    
stidx = 66000
edidx = 69000
yD = yD[stidx:edidx]
cyD = cyD[stidx:edidx]
maxyD = max(yD)
minyD = min(yD)

#(idxmyD,idxMyD)

yO = yO[stidx:edidx]
cyO = cyO[stidx:edidx]
maxyO = max(yO)
minyO = min(yO)



fig = plt.figure(figsize=(20,14),constrained_layout = True)
#ax = matplotlib.axes._axes.Axes[6]
gs = fig.add_gridspec(8,8)
f3_ax1 = fig.add_subplot(gs[2:4,0:2])
f3_ax2 = fig.add_subplot(gs[2:3,2:4])
f3_ax3 = fig.add_subplot(gs[3:4,2:4])

f3_ax4 = fig.add_subplot(gs[0:2,0:4])
f3_ax5 = fig.add_subplot(gs[0:2,4:8])
f3_ax6 = fig.add_subplot(gs[2:4,4:8])
f3_ax7 = fig.add_subplot(gs[4:6,4:8])
f3_ax8 = fig.add_subplot(gs[4:6,0:4])
f3_ax9 = fig.add_subplot(gs[6:8,0:4])
f3_ax10 = fig.add_subplot(gs[6:8,4:8])




ax = fig.axes
ax[0] = f3_ax4
ax[1] = f3_ax5
ax[2] = f3_ax6
ax[3] = f3_ax7
ax[4] = f3_ax8
ax[5] = f3_ax9
ax[6] = f3_ax10



f3_ax2.plot(t[stidx:edidx],yD[0:edidx-stidx+1],color = 'xkcd:black',linestyle = 'dashed', lw = 0.7)
f3_ax2.plot(t[stidx:edidx],cyD[0:edidx-stidx+1],color = 'xkcd:black',linestyle = 'solid', lw = 0.6)
f3_ax3.plot(t[stidx:edidx],yO[0:edidx-stidx+1],color = 'xkcd:blue',linestyle = 'dashed', lw = 0.7)
f3_ax3.plot(t[stidx:edidx],cyO[0:edidx-stidx+1],color = 'xkcd:blue',linestyle = 'solid', lw = 0.6)



f3_ax1.plot(t[stidx:edidx],yD[0:edidx-stidx+1],color = 'xkcd:black',linestyle = 'solid',linewidth = 0.7)
f3_ax1.plot(t[stidx:edidx],yO[0:edidx-stidx+1],color = 'xkcd:blue',linestyle = 'solid',linewidth = 0.7)
f3_ax1.set_xlim([662,688])
f3_ax1.set_ylim([-0.5,0.5])


xmajorLocator = MultipleLocator(5.0)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(1.0)
f3_ax1.xaxis.set_minor_locator(xminorLocator)
f3_ax1.xaxis.set_major_locator(xmajorLocator)
f3_ax1.set_xlabel('$time$',fontsize=10,color='xkcd:black')
f3_ax1.set_ylabel('$Y$',fontsize=10,color='xkcd:black')

f3_ax1.tick_params(which = 'both',left=1,top=1,right=1,
                  bottom=1,labelleft=1,labeltop=0,
                  labelright=0,labelbottom=1,width=0.5,direction="in")
f3_ax1.text(-0.06, 0.93, '$(i)$',
        verticalalignment='top', horizontalalignment='right',
        transform=f3_ax1.transAxes,
        color='black', fontsize=12)


f3_ax2.xaxis.set_minor_locator(xminorLocator)
f3_ax2.xaxis.set_major_locator(xmajorLocator)
f3_ax2.set_ylabel('$C_y$',fontsize = 10, color='xkcd:black')
f3_ax2.yaxis.set_major_locator(MultipleLocator(1.0))
f3_ax2.yaxis.set_minor_locator(MultipleLocator(0.25))
ax2 = f3_ax2.twinx()
ax2.set_ylabel('$y$')
f3_ax2.set_xlabel('$time$',fontsize = 10, color='xkcd:black')
f3_ax2.tick_params(which = 'both',left=1,top=1,right=1,
                  bottom=1,labelleft=1,labeltop=0,
                  labelright=0,labelbottom=1,width=0.5,direction="in")
ax2.tick_params(which = 'both',left=0,top=0,right=1,
                  bottom=0,labelleft=0,labeltop=0,
                  labelright=0,labelbottom=0,width=0.5,direction="in")
f3_ax2.set_xlim([662,688])
f3_ax2.set_ylim([-2.5,2.5])
f3_ax2.text(-0.06, 0.93, '$(ii)$',
        verticalalignment='top', horizontalalignment='right',
        transform=f3_ax2.transAxes,
        color='black', fontsize=12)



f3_ax3.xaxis.set_minor_locator(xminorLocator)
f3_ax3.xaxis.set_major_locator(xmajorLocator)
f3_ax3.set_ylabel('$C_y$',fontsize = 10, color='xkcd:black')
f3_ax3.yaxis.set_major_locator(MultipleLocator(1.0))
f3_ax3.yaxis.set_minor_locator(MultipleLocator(0.2))
ax2 = f3_ax3.twinx()
ax2.set_ylabel('$y$')
f3_ax3.set_xlabel('$time$',fontsize = 10, color='xkcd:black')
f3_ax3.tick_params(which = 'both',left=1,top=1,right=1,
                  bottom=1,labelleft=1,labeltop=0,
                  labelright=0,labelbottom=1,width=0.5,direction="in")
ax2.tick_params(which = 'both',left=0,top=0,right=0,
                  bottom=0,labelleft=0,labeltop=0,
                  labelright=0,labelbottom=0,width=0.5,direction="in")
f3_ax3.set_xlim([662,688])
f3_ax3.set_ylim([-2.5,1.7])
f3_ax3.text(-0.06, 0.93, '$(iii)$',
        verticalalignment='top', horizontalalignment='right',
        transform=f3_ax3.transAxes,
        color='black', fontsize=12)




e = mpatches.Ellipse((669.5,maxyD),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(669.5,maxyD,'a')
e = mpatches.Ellipse((670.5,minyO),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(670.5,minyO-0.05,'b')
e = mpatches.Ellipse((671.5,minyD),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(671.5,minyD,'c')
e = mpatches.Ellipse((672.5,maxyO),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(672.5,maxyO,'d')
e = mpatches.Ellipse((673.5,maxyD),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(673.5,maxyD,'e')
e = mpatches.Ellipse((674.5,maxyD),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(674.5,maxyD,'f')
e = mpatches.Ellipse((675.5,maxyD),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(675.5,maxyD,'g')


labels = ['a','b','c','d','e','f','g']


folder_path0 = 'I:\\NIHAR\\Tandem_D_O\\D_O\\D_O_3.5\\plt'

pltfilelist = ['q0066950.plt','q0067050.plt','q0067150.plt','q0067250.plt','q0067350.plt','q0067450.plt','q0067550.plt']
for i in range(len(pltfilelist)):
    A1 = os.path.join(folder_path0,pltfilelist[i])
    
    xc,yc,qqf = rvd.read_binary(A1)
    X,Y = np.meshgrid(xc,yc)
    qqf[:,:,0,4] = np.minimum(qqf[:,:,0,4],2.0)
    qqf[:,:,0,4] = np.maximum(qqf[:,:,0,4],-2.0)
    
    mask = (qqf[:,:,0,5] > 0.5)
    vort = ma.masked_array(qqf[:,:,0,4],mask = mask)
    
    ax[i].contourf(X,Y,vort,10,cmap=plt.cm.bwr)
    xlim = [4.5,22.5]
    ylim = [7.5,12.5]
    ax[i].set_xlim(xlim)
    ax[i].set_ylim(ylim)
    ax[i].set_aspect(1)
    ax[i].set_adjustable("box")
    ax[i].set_xticks([])
    ax[i].set_yticks([])
    
    ax[i].plot([4.7,4.7],[10.0+minyD,10.0+maxyD],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    ax[i].plot([4.65,4.75],[10.0+maxyD,10.0+maxyD],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    ax[i].plot([4.65,4.75],[10.0+minyD,10.0+minyD],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    
    ax[i].plot([9.7,9.7],[10.0+minyO,10.0+maxyO],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    ax[i].plot([9.65,9.75],[10.0+minyO,10.0+minyO],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    ax[i].plot([9.65,9.75],[10.0+maxyO,10.0+maxyO],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    
    ax[i].text(0.02, 0.93, labels[i],
        verticalalignment='top', horizontalalignment='right',
        transform=ax[i].transAxes,
        color='black', fontsize=12)
    
fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Vort_3.5.pdf', dpi=300,format='pdf',bbox_inches='tight', pad_inches=0.1)





