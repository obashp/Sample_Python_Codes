import os
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import gridspec
import Atul_myplt as myplt
import matplotlib.lines as mlines

plt.close('all')


majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1)

majorLocator1 = MultipleLocator(1)
majorFormatter1 = FormatStrFormatter('%d')
minorLocator1 = MultipleLocator(0.2)


rowNum = 2  # num of rows of the subplots
colNum = 1  # num of colums of the subplots

figW, figH = myplt.figsize_mm(figSWmm=300, figAR=0.4,
                              nRows=rowNum, nCols=colNum)

#fig, axs = plt.subplots(figsize=(figW, figH),
#                        nrows=rowNum, ncols=colNum,
#                        sharex=False,squeeze=False)  # sharing the x-label
#


fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=rowNum, ncols=colNum,
                        sharex=False,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=5,
                    wspace=2,hspace=4)


folder_path0 = 'F:\\Nihar\\VIV_tests\\D_O'


result = 'post'

Urfile = os.path.join(folder_path0,'Urlist.dat')

with open(Urfile, 'r') as f:
    Ur = np.loadtxt(f)
    
Ur.sort()

NI = len(Ur)
NJ = 301  
NUMVAR = 2  
xc = np.zeros((NI,1),dtype=float)
yc = np.zeros((NJ,1),dtype=float)
qq = np.zeros((NJ,NI,NUMVAR),dtype=float)

#St0 = 0.128174 # Re = 50

St0 = 0.168067#St for D_section, W=8D at Re=100
print(St0)
#St0 = 0.2838 #St for D_section, W=8D at Re=100

#Ur = [2., 3., 4., 5., 6., 7., 8., 9. ,10.,11.,12.,13.,15.,18.,20.,21.,22.,25.,26.,28.,29.,30.,50.,100.]
#Ur  = [2.,3.,4.,5.,6.,8.,10.,12.,14.,16.,18.,20.]
#Ur = [2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.]

St1 = np.zeros(NI)
fN = np.zeros(NI)


for i in range(NI):
        #St1[i] = St0*Ur[i]*np.sqrt(1.5) used for added mass effect
    St1[i] = St0*Ur[i]
    fN[i] = 1.0
#plt.cla()
    
ax1 = axs[0,0]
ax2 = axs[1,0]
#*************************************************************************************************
filename = os.path.join(folder_path0,'PSD_cly_body2.plt')



with open(filename, 'r') as f: 
    data1 = np.loadtxt(f,skiprows=2)
 

x = data1[:,0]
y = data1[:,1]

#print(x)

xll = x.min();  xul = x.max();  yll = y.min();  yul = y.max()

xc = Ur #np.linspace(xll, xul, NI)
yc = np.linspace(yll, yul, NJ) 

    
X,Y = np.meshgrid(xc,yc)

count = 0
for i in range(NI):
    for j in range(NJ):
                qq[j,i,0] = data1[count,3]
                qq[j,i,1] = data1[count,4]
                count = count + 1


qq[:,:,0] = np.minimum(qq[:,:,0],3)
qq[:,:,0] = np.maximum(qq[:,:,0],0)


def safe_ln(x, minval=0.001):
    return np.log(x.clip(min=minval))


fq = qq[:,:,0]

fq = safe_ln(fq)

ax1.contourf(X,Y,fq,50,cmap=plt.cm.afmhot_r) 


ax1.plot(Ur,St1, linestyle='--', label = '', color = 'b')
ax1.plot(Ur,fN, linestyle='--', label = '', color = 'b')

#ax.axis('scaled')    
ax1.set_xlim([1,20])
ax1.set_ylim([0,6])
 
fonts = 13

ax1.set_xticklabels([])
#ax.set_xlabel('$U_{R}$', fontsize=fonts, labelpad=0.1)
ax1.set_ylabel('$f_{Cy}$', fontsize=fonts, color = 'black')
    
ax1.xaxis.set_minor_locator(minorLocator)
ax1.xaxis.set_major_locator(majorLocator)
ax1.yaxis.set_minor_locator(minorLocator1)
ax1.yaxis.set_major_locator(majorLocator1)


ax1.text(-0.05, 0.93, '(a)',
        verticalalignment='top', horizontalalignment='right',
        transform=ax1.transAxes,
        color='black', fontsize=fonts)

#ax1.text(0.69, 0.35, '$f_{v_D}$ =  0.181076',
#        verticalalignment='top', horizontalalignment='right',
#        transform=ax1.transAxes,
#        color='black', fontsize=fonts+3)

#ax.text(0.15, 0.9, '$Ri$ = 0',
#        verticalalignment='top', horizontalalignment='right',
#        transform=ax.transAxes,
#        color='black', fontsize=fonts)
    
ax1.tick_params(which='both',left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=0, width=0.5, direction="in")
ax1.tick_params(which='major', left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=0,length=5, direction="in")
ax1.tick_params(which='minor', left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=0,length=3,direction="in")

qq[:,:,1] = np.minimum(qq[:,:,1],3)
qq[:,:,1] = np.maximum(qq[:,:,1],0)


fq = qq[:,:,1]

fq = safe_ln(fq)


ax = ax2.contourf(X,Y,fq,50,cmap=plt.cm.afmhot_r) 
ax2.plot(Ur,St1, linestyle='--', label = '$Ri = 0$', color = 'b')
ax2.plot(Ur,fN, linestyle='--', label = '$Ri = 0$', color = 'b')
#ax.axis('scaled')    
ax2.set_xlim([1,20])
ax2.set_ylim([0,6])

#ax2.set_xlabel('$U_{R}$', fontsize=fonts, labelpad=0.1)
ax2.set_ylabel('$f_{Cv}$', fontsize=fonts, color = 'black')
   
    
ax2.xaxis.set_minor_locator(minorLocator)
ax2.xaxis.set_major_locator(majorLocator)
ax2.yaxis.set_minor_locator(minorLocator1)
ax2.yaxis.set_major_locator(majorLocator1)

ax2.text(-.05, 0.93, '(b)',
        verticalalignment='top', horizontalalignment='right',
        transform=ax2.transAxes,
        color='black', fontsize=fonts) 

#ax1.text(0.15, 0.9, '$Ri$ = 0',
#        verticalalignment='top', horizontalalignment='right',
#        transform=ax1.transAxes,
#        color='black', fontsize=fonts)
    
ax2.tick_params(which='both',left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=1, width=0.5, direction="in")
ax2.tick_params(which='major', left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=1,length=5, direction="in")
ax2.tick_params(which='minor', left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=0,length=3,direction="in")

#ax2.set_xlabel('$U_R$', fontsize=12)


fig.colorbar(ax, orientation='horizontal')
fig.tight_layout()
fig.savefig('PSD_O.jpg', dpi=150,format='jpg',bbox_inches='tight', pad_inches=0.1)
