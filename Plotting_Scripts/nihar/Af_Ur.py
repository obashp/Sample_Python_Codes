import os
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import gridspec
from matplotlib.patches import Ellipse
import Atul_myplt as myplt
import matplotlib.lines as mlines

plt.close('all')


def safe_ln(x, minval=0.001):
    return np.log(x.clip(min=minval))


majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1)

majorLocator1 = MultipleLocator(1)
majorFormatter1 = FormatStrFormatter('%d')
minorLocator1 = MultipleLocator(0.2)


rowNum = 4  # num of rows of the subplots
colNum = 1  # num of colums of the subplots

figW, figH = myplt.figsize_mm(figSWmm=290, figAR=0.33,
                              nRows=rowNum, nCols=colNum)


fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=rowNum, ncols=colNum,
                        sharex=True,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=10, right=10, top=5, bottom=5,
                    wspace=2,hspace=2)


folder_path0 = 'F:\\Nihar\\VIV_tests\\D_O'


result = 'post'

Urfile = os.path.join(folder_path0,'Urlist.dat')

with open(Urfile, 'r') as f:
    Ur = np.loadtxt(f)
    
Ur.sort()

NI = len(Ur)
NJ = 301  
NUMVAR = 3  
xc = np.zeros((NI,1),dtype=float)
yc = np.zeros((NJ,1),dtype=float)
qq = np.zeros((NJ,NI,NUMVAR),dtype=float)



St0 = 0.181736#St for D_section, W=8D at Re=100
print(St0)

St1 = np.zeros(NI)
fN = np.zeros(NI)


labels = ['$f_{y}$','$f_{Cy}$','$f_{Cv}$','(a)','(b)','(c)','$f_{Cx}$','(d)']
fonts = 15


for i in range(NI):
    St1[i] = St0*Ur[i]
    fN[i] = 1.0
    
filename = os.path.join(folder_path0,'PSD_cly_body1.plt')

with open(filename, 'r') as f: 
    data1 = np.loadtxt(f,skiprows=2)
 
x = data1[:,0]
y = data1[:,1]

yll = y.min();  yul = y.max()
xc = Ur 
yc = np.linspace(yll, yul, NJ) 
X,Y = np.meshgrid(xc,yc)

levels = [-3,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]

count = 0
for i in range(NI):
    for j in range(NJ):
        qq[j,i,0] = data1[count,2]
        qq[j,i,1] = data1[count,3]
        qq[j,i,2] = data1[count,4]
        count = count+1

           
for t in range(3):
    qq[:,:,t] = np.minimum(qq[:,:,t],3)
    qq[:,:,t] = np.maximum(qq[:,:,t],0)
    fq = qq[:,:,t]
    fq = safe_ln(fq)
    
    ax = axs[t,0]
    ax1 = ax.contourf(X,Y,fq,levels,cmap=plt.cm.afmhot_r)
   
    if t==0:
        ax.text(0.59, 0.45, '$f_{v_D}$ =  ' + str(St0),
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=fonts+3)
        
    
        
    ax.plot(Ur,St1, linestyle='--', label = '', color = 'k')
    ax.plot(Ur,fN, linestyle='--', label = '', color = 'grey')
    ax.plot(Ur,3*fN, linestyle='--', label = '', color = 'grey')
    ax.plot(Ur,5*fN, linestyle='--', label = '', color = 'grey')
    ax.set_xlim([1,30])
    ax.set_ylim([0,6])
    
    

    ax.set_ylabel(labels[t], fontsize=fonts, color = 'black')
    
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator1)
    ax.yaxis.set_major_locator(majorLocator1)
    
    ax.text(-0.05, 0.93, labels[t+3],
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=fonts)
    
    ax.tick_params(which='both',left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=0, width=0.5, direction="in")
    ax.tick_params(which='major', left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=(1 if t==2 else 0),length=5, direction="in")
    ax.tick_params(which='minor', left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=0,length=3,direction="in")
    
 
filename = os.path.join(folder_path0,'PSD_cdx_body1.plt')

with open(filename, 'r') as f: 
    data1 = np.loadtxt(f,skiprows=2)
 
x = data1[:,0]
y = data1[:,1]

yll = y.min();  yul = y.max()
xc = Ur 
yc = np.linspace(yll, yul, NJ) 
X,Y = np.meshgrid(xc,yc)    
count = 0
for i in range(NI):
    for j in range(NJ):
        qq[j,i,0] = data1[count,2]
        count = count+1
    
qq[:,:,0] = np.minimum(qq[:,:,0],3)
qq[:,:,0] = np.maximum(qq[:,:,0],0)
fq = qq[:,:,0]
fq = safe_ln(fq)
    
ax = axs[3,0]
ax1 = ax.contourf(X,Y,fq,levels,cmap=plt.cm.afmhot_r)
   
      
    
        
ax.plot(Ur,fN, linestyle='--', label = '', color = 'grey')
ax.plot(Ur,3*fN, linestyle='--', label = '', color = 'grey')
ax.plot(Ur,5*fN, linestyle='--', label = '', color = 'grey')
ax.set_xlim([1,30])
ax.set_ylim([0,6])
    
ax.set_ylabel(labels[6], fontsize=fonts, color = 'black')
    
ax.xaxis.set_minor_locator(minorLocator)
ax.xaxis.set_major_locator(majorLocator)
ax.yaxis.set_minor_locator(minorLocator1)
ax.yaxis.set_major_locator(majorLocator1)
    
ax.text(-0.05, 0.93, labels[7],
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=fonts)
    
ax.tick_params(which='both',left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=0, width=0.5, direction="in")
ax.tick_params(which='major', left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=(1 if t==2 else 0),length=5, direction="in")
ax.tick_params(which='minor', left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=0,length=3,direction="in")
   
ax.set_xlabel('$U_R$', fontsize=fonts)
fig.colorbar(ax1, ax = axs[3,0],orientation='horizontal',fraction=0.1, aspect=40,shrink=0.8)
fig.savefig('C:\\Users\\Laxman\\Documents\\POF1\PSD_D.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)





