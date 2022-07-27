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
import SciUtils as Sut
from matplotlib.patches import Ellipse,Arc


def marker_plot(plotloc,filename):
    Npoint = []
    lineidx = []
    k = 0
    with open(filename,'r') as f:
        while(1):
            line = f.readline()
            if line == '':
                break
            line = line.split()
            # print(line)
            if(line[0] == 'ZONE'):
                Npoint.append(int(int(line[2])/3))
                lineidx.append(k)
        
            k = k+1
    
    body = 0
    while body < len(Npoint):
        with open(markerfile,'r') as f:
            data = np.loadtxt(f,dtype=float,skiprows=lineidx[body]+1,max_rows=int(Npoint[body]))
            x = data[:,0]
            y = data[:,1]
            plotloc.plot(x,y,c='xkcd:black',lw=0.8)
        body = body+1
        
    
   


folder_path = 'I:\\NIHAR\\Tandem_D_O\\D_O\\D_O_3.5'

fig = plt.figure(figsize=(20,14),constrained_layout = True)
# gs = fig.add_gridspec(10,4)
# f3_ax1 = fig.add_subplot(gs[0:2,0:2])
# f3_ax2 = fig.add_subplot(gs[0:1,2:4])
# f3_ax3 = fig.add_subplot(gs[1:2,2:4])
# f3_ax4 = fig.add_subplot(gs[2:4,0:2])
# f3_ax5 = fig.add_subplot(gs[2:4,2:4])
# f3_ax6 = fig.add_subplot(gs[4:6,0:2])
# f3_ax7 = fig.add_subplot(gs[4:6,2:4])
# f3_ax8 = fig.add_subplot(gs[6:8,0:2])
# f3_ax9 = fig.add_subplot(gs[6:8,2:4])
# f3_ax10 = fig.add_subplot(gs[8:10,0:2])
# f3_ax11 = fig.add_subplot(gs[8:10,2:4])

gs = fig.add_gridspec(10,4)
f3_ax1 = fig.add_subplot(gs[0:2,0:2])
f3_ax2 = fig.add_subplot(gs[0:1,2:4])
f3_ax3 = fig.add_subplot(gs[1:2,2:4])
f3_ax4 = fig.add_subplot(gs[2:4,0:2])
f3_ax5 = fig.add_subplot(gs[2:4,2:4])
f3_ax6 = fig.add_subplot(gs[4:6,0:2])
f3_ax7 = fig.add_subplot(gs[4:6,2:4])
# f3_ax8 = fig.add_subplot(gs[6:8,0:2])
# f3_ax9 = fig.add_subplot(gs[6:8,2:4])
# f3_ax10 = fig.add_subplot(gs[8:10,0:2])
# f3_ax11 = fig.add_subplot(gs[8:10,2:4])

ax = fig.axes
ax[0] = f3_ax4
ax[1] = f3_ax5
ax[2] = f3_ax6
ax[3] = f3_ax7
# ax[4] = f3_ax8
# ax[5] = f3_ax9
# ax[6] = f3_ax10
# ax[7] = f3_ax11




tlim,ylim,cylimD,cylimO,AyD,AyO = Sut.plot_ty(f3_ax1,f3_ax2,f3_ax3,folder_path)

yfval = round(0.2*(ylim[1]-ylim[0]),1)
ut.plotparams(f3_ax1,'$y$','$time$',tlim,ylim,True,True,yfval,True,True,5.0)
f3_ax1.legend(['$y_O$','$y_D$'],edgecolor='1.0',loc='upper right')
yfval = round(0.2*(cylimD[1]-cylimD[0]),1)
ut.plotparams(f3_ax2,'$C_y$','$time$',tlim,cylimD,True,True,yfval,True,True,5.0)
ax1 = f3_ax2.twinx()
ax1.set_ylabel('$y$')
f3_ax2.legend(['$Cy$','$y$'],edgecolor='1.0',loc='upper right')
yfval = round(0.2*(cylimO[1]-cylimO[0]),1)
ut.plotparams(f3_ax3,'$C_y$','$time$',tlim,cylimO,True,True,yfval,True,True,5.0)
ax1 = f3_ax3.twinx()
ax1.set_ylabel('$y$')
f3_ax3.legend(['$Cy$','$y$'],edgecolor='1.0',loc='upper right')

folder_path1 = folder_path + '\\plt\\'
folder_path2 = folder_path + '\\marker\\'

# UR 3.5
pltfilelist = ['q0172250.plt','q0172450.plt','q0172700.plt','q0172850.plt']
markerfilelist = ['marker.0172250.dat','marker.0172450.dat','marker.0172700.dat','marker.0172850.dat']
# pltfilelist = ['q0172000.plt','q0172100.plt','q0172200.plt','q0172300.plt','q0172400.plt','q0172500.plt','q0172600.plt','q0172650.plt']
# markerfilelist = ['marker.0172000.dat','marker.0172100.dat','marker.0172200.dat','marker.0172300.dat','marker.0172400.dat','marker.0172500.dat','marker.0172600.dat','marker.0172650.dat']

# UR 4
# pltfilelist = ['q0080300.plt','q0080400.plt','q0080500.plt','q0080650.plt','q0080700.plt','q0080850.plt']
# markerfilelist = ['marker.0080300.dat','marker.0080400.dat','marker.0080500.dat','marker.0080650.dat','marker.0080700.dat','marker.0080850.dat']

# UR 4.5
# pltfilelist = ['q0100250.plt','q0100350.plt','q0100400.plt','q0100500.plt','q0100600.plt','q0100700.plt']
# markerfilelist = ['marker.0100250.dat','marker.0100350.dat','marker.0100400.dat','marker.0100500.dat','marker.0100600.dat','marker.0100700.dat']

# UR 6.5
# pltfilelist = ['q0063650.plt','q0063700.plt','q0063800.plt','q0064050.plt','q0064150.plt','q0064300.plt']
# markerfilelist = ['marker.0063650.dat','marker.0063700.dat','marker.0063800.dat','marker.0064050.dat','marker.0064150.dat','marker.0064300.dat']

# UR 8.5
# pltfilelist = ['q0096050.plt','q0096250.plt','q0096700.plt','q0097000.plt','q0097200.plt','q0097400.plt','q0097550.plt','q0097850.plt']
# markerfilelist = ['marker.0096050.dat','marker.0096250.dat','marker.0096700.dat','marker.0097000.dat','marker.0097200.dat','marker.0097400.dat','marker.0097550.dat','marker.0097850.dat']

# UR 12.5
# pltfilelist = ['q0080100.plt','q0080250.plt','q0080300.plt','q0080450.plt','q0080550.plt','q0080650.plt']
# markerfilelist = ['marker.0080100.dat','marker.0080250.dat','marker.0080300.dat','marker.0080450.dat','marker.0080550.dat','marker.0080650.dat']


labels = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)']
for i in range(len(pltfilelist)):
    A1 = os.path.join(folder_path1,pltfilelist[i])
    markerfile = os.path.join(folder_path2,markerfilelist[i])
    
    xc,yc,qqf = rvd.read_binary(A1)
    X,Y = np.meshgrid(xc,yc)
    qqf[:,:,0,4] = np.minimum(qqf[:,:,0,4],2.0)
    qqf[:,:,0,4] = np.maximum(qqf[:,:,0,4],-2.0)
    
    
    mask = (qqf[:,:,0,5] > 0.5)
    vort = ma.masked_array(qqf[:,:,0,4],mask = mask)
    
    ax[i].contourf(X,Y,vort,10,cmap=plt.cm.bwr)
    xlim = [4.5,23.0]
    ylim = [7.5,12.5]
    # ylim = [6.5,13.5]
    # ylim = [5.0,15.0]
    ax[i].set_xlim(xlim)
    ax[i].set_ylim(ylim)
    ax[i].set_aspect(1)
    ax[i].set_adjustable("box")
    ax[i].set_xticks([])
    ax[i].set_yticks([])
    
    marker_plot(ax[i],markerfile)
    
    ax[i].plot([4.7,4.7],[10.0+AyD[0],10.0+AyD[1]],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    ax[i].plot([4.65,4.75],[10.0+AyD[1],10.0+AyD[1]],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    ax[i].plot([4.65,4.75],[10.0+AyD[0],10.0+AyD[0]],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    
    ax[i].plot([9.7,9.7],[10.0+AyO[0],10.0+AyO[1]],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    ax[i].plot([9.65,9.75],[10.0+AyO[0],10.0+AyO[0]],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    ax[i].plot([9.65,9.75],[10.0+AyO[1],10.0+AyO[1]],color='xkcd:gray',linestyle='solid',linewidth = 1.5)
    
    ax[i].text(-0.05, 0.93, labels[i],
        verticalalignment='top', horizontalalignment='right',
        transform=ax[i].transAxes,
        color='black', fontsize=12)
    
    

c = Ellipse(xy=(7.7, 10.0), width=4.0, height=4,ls = 'dashed', 
            edgecolor='k', fc='None', lw=1)
ax[0].add_patch(c)
c = Ellipse(xy=(18.9, 9.0), width=5.5, height=3.0,angle = -60,ls = 'dashed', 
            edgecolor='k', fc='None', lw=1)
ax[0].add_patch(c)

c = Ellipse(xy=(7.6, 10.0), width=4.0, height=2.3,angle = -150,ls = 'dashed', 
            edgecolor='k', fc='None', lw=1)
ax[1].add_patch(c)
c = Ellipse(xy=(12.8, 10.3), width=4.5, height=2.0,angle = -150,ls = 'dashed', 
            edgecolor='k', fc='None', lw=1)
ax[1].add_patch(c)

c = Ellipse(xy=(7, 10.0), width=3.0, height=3.0,angle = 105,ls = 'dashed', 
            edgecolor='k', fc='None', lw=1)
ax[2].add_patch(c)
c = Ellipse(xy=(14.3, 10.3), width=5.0, height=2.3,angle = 0,ls = 'dashed', 
            edgecolor='k', fc='None', lw=1)
ax[2].add_patch(c)
c = Ellipse(xy=(11.7, 8.7), width=2.0, height=2.0,angle = 0,ls = 'dashed', 
            edgecolor='k', fc='None', lw=1)
ax[2].add_patch(c)

e = mpatches.Ellipse((1722.5,0.35),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
# f3_ax2.add_artist(e)
# f3_ax3.add_artist(e)
f3_ax1.text(1722.5,0.24,'a')
e = mpatches.Ellipse((1724.5,0.35),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(1724.5,0.21,'b')
e = mpatches.Ellipse((1727.0,0.35),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(1727.0,0.21,'c')
e = mpatches.Ellipse((1728.5,0.32),0.4,0.04,edgecolor='k',fc='None', lw=0.5)
f3_ax1.add_artist(e)
f3_ax1.text(1728.5,0.32,'d')


fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Vort_3.5.pdf', dpi=300,format='pdf',bbox_inches='tight', pad_inches=0.1)


