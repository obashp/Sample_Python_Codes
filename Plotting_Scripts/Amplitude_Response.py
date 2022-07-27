import os
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import Atul_myplt as myplt
import Utils as Ut
from matplotlib.patches import Ellipse

def plotDO(plotloc,aspect,parameters):
    r = parameters[2]
    xc = parameters[0]
    yc = parameters[1]
    
    pi = 4*math.atan(1.0)
    npts = 181
    dtheta = pi/(npts-1)
    i = 0
    
    coord = np.zeros((npts,2),dtype=float)
    
    while i < npts:
        coord[i][0] = xc + r*math.cos(-pi/2 + i*dtheta)/aspect
        coord[i][1] = yc + r*math.sin(-pi/2 + i*dtheta)
        i = i+1
    x,y = zip(*coord)
    plotloc.plot(x,y,color=parameters[3])
    plotloc.plot([x[0],x[0]],[y[0],y[-1]],color=parameters[3])
     
    i = 0
    coord = np.zeros((361,2),dtype=float)
    dtheta = pi/180
    while i < 361:
        coord[i][0] = xc + 1.0 + r*math.cos(i*dtheta)/aspect
        coord[i][1] = yc + r*math.sin(i*dtheta)
        i = i+1
    #print(coord)
    x,y = zip(*coord)
    plotloc.plot(x,y,color=parameters[4])


folder_path = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\D_O\\'

body = '2'
dire = 'y'
St = 0.1680


A = '_AC_Ur.dat'
PSD = 'PSD_body_'
Ampfile = os.path.join(folder_path,body + dire +A)
PSDfile = os.path.join(folder_path,PSD + body + dire +'.plt')

Ur = Ut.GetUrspan(folder_path)
nUr = len(Ur)
#print(nUr)



figW, figH = myplt.figsize_mm(figSWmm=200, figAR=0.19,
                              nRows=6, nCols=1)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=6, ncols=1,
                        sharex=True,squeeze=False)

myplt.margins_adjust(fig, figW, figH, left=10, right=10, top=5, bottom=5,
                   wspace=0,hspace=2)

parameters = [False,True,0.25]
ax = axs[0,0]
aspect = Ut.Amplitudeplot(ax,Ampfile,nUr,'kv-', '(a)',False,True,0.4)
print(aspect)

#circularcylinder = [0.5,2,1.1,0.15,'black','red']
#Dsection = [0.9,2,2.1,0.25,'red','black']

# plotDO(ax,0.2,[2,0.3,0.05,'xkcd:black','xkcd:red'])
plotDO(ax,0.5,[2,1.1,0.15,'xkcd:black','xkcd:red'])
# plotDO(ax,0.9,[2,2.1,0.25,'xkcd:red','xkcd:black'])
ax = axs[1,0]
ylabel = '$f_' + dire +'$'
Ut.plotStfN(ax,True,St,True,Ur)
ax1 = Ut.PSDplot(ax, PSDfile, Ur, St, 2, ylabel,'(b)')
#ax.set_aspect(1)
ax = axs[2,0]
ylabel = '$f_{C' + dire+'}$'
Ut.plotStfN(ax,True,St,False,Ur)
# Ut.plotStfN(ax,True,0.168067,False,Ur)
ax1 = Ut.PSDplot(ax, PSDfile, Ur, St, 3, ylabel,'(c)')
ax.arrow(10.0,2.5,2.5,0.315,head_width=0.1)
ax.arrow(10.0,1.6,2.5,-0.315,head_width=0.1)

ax = axs[3,0]
ylabel = '$f_{C' + dire+'v}$'
Ut.plotStfN(ax,True,St,False,Ur)
# Ut.plotStfN(ax,True,0.168067,False,Ur)
ax1 = Ut.PSDplot(ax, PSDfile, Ur, St, 4, ylabel,'(d)')
ax.arrow(10.0,2.6,2.5,0.315,head_width=0.1)
ax.arrow(10.0,1.61,2.5,-0.315,head_width=0.1)
#ax.set_aspect(1)
# ax = axs[3,0]
# dire = 'x'
# ylabel = '$f_{C' + dire +'}$'
# PSDfile = os.path.join(folder_path,PSD + body + dire +'.plt')
# ax1 = Ut.PSDplot(ax, PSDfile, Ur, St, 2, ylabel, '(d)', True)

# fig.subplots_adjust(bottom = 0.1)
# cbar_ax = fig.add_axes([0.25, 0.03, 0.5, 0.02])
# fig.colorbar(ax1, cax = cbar_ax, orientation='horizontal')

A = 'phi_Ur.dat'
Ampfile = os.path.join(folder_path,body + A)
ax = axs[4,0]
ylabel = '$\u03A6_{C_' + dire+'}$'
ax1 = Ut.PSDplot(ax,Ampfile,Ur,St,2,ylabel,'(e)',False,False)

ax = axs[5,0]
ylabel = '$\u03A6_{C_{v' + dire+'}}$'
ax1 = Ut.PSDplot(ax,Ampfile,Ur,St,3,ylabel,'(f)',True,False)

#fig.colorbar(ax1, ax = axs[3,0],orientation='horizontal',fraction=0.1, aspect=40,shrink=0.8)
fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\Amp_yO.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)