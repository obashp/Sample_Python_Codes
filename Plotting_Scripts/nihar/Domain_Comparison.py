import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Atul_myplt as myplt
import Utils as Ut
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes,mark_inset,InsetPosition)


folder_path = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\Grid_Domain\\Ur_8\\'

colors = ['xkcd:black','xkcd:blue','xkcd:red','xkcd:green','xkcd:gray']
ls = ['-','-.','--',':','-']
timeaddD = [0.0,0.32,-0.25,-0.2,0.25]#np.zeros(5,dtype = float)
timeaddO = [0.0,0.3,-0.27,-0.35,0.0]#np.zeros(5,dtype = float)

#folder_list = ['35x20_5D','35x20_5D','35x20_5D','35x20_5D']
# y = [x[0] for x in os.walk(folder_path)]
# print(y)

foldernames = [f.name for f in os.scandir(folder_path) if f.is_dir()]
subfolders = [f.path for f in os.scandir(folder_path) if f.is_dir()]

figW, figH = myplt.figsize_mm(figSWmm=250, figAR=0.5,
                              nRows=2, nCols=1)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=2, ncols=1,
                        sharex=True,squeeze=False)

myplt.margins_adjust(fig, figW, figH, left=10, right=10, top=5, bottom=5,
                   wspace=0,hspace=2)

ax = axs[0,0]
bx = axs[1,0]
axins = plt.axes([0,0,1,1])
bxins = plt.axes([5,0,1,1])
ap = InsetPosition(ax,[0.45,0.6,0.3,0.3])
bp = InsetPosition(bx,[0.45,0.65,0.3,0.3])
axins.set_axes_locator(ap)
bxins.set_axes_locator(bp)
mark_inset(ax,axins,loc1=3,loc2=4,fc="none",ec='0.5')
mark_inset(bx,bxins,loc1=3,loc2=4,fc="none",ec='0.5')
Ut.plotparams(bxins, '','',[573.5,575], [0.84,1.05],True,True,0.04,True,True,0.5)
Ut.plotparams(axins, '','',[574.0,575.2], [1.38,1.5],True,True,0.025,True,True,0.5)


i = 0
while i < len(subfolders):
    Afile = os.path.join(subfolders[i],'data_'+foldernames[i]+'.dat')
    print(foldernames[i])
    with open(Afile,'r') as f:
        data = np.loadtxt(f,skiprows=1)
        time = data[:,0]
        y1 = data[:,1]
        y2 = data[:,2]
    
    time1 = time + timeaddD[i]
    sti = np.searchsorted(time,572)
    edi = np.searchsorted(time,576)
    ax.plot(time1,y1,colors[i],linestyle=ls[i])
    axins.plot(time1[sti:edi],y1[sti:edi],colors[i],linestyle=ls[i])
    
    time = time + timeaddO[i]
    sti = np.searchsorted(time,573.5)
    edi = np.searchsorted(time,575.2)
    bx.plot(time,y2,colors[i],linestyle=ls[i])
    bxins.plot(time[sti:edi],y2[sti:edi],colors[i],linestyle=ls[i])
    i = i+1

ax.text(-0.02, 0.93, '(a)',
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=12)

bx.text(-0.02, 0.93, '(b)',
        verticalalignment='top', horizontalalignment='right',
        transform=bx.transAxes,
        color='black', fontsize=12)
    

ax.legend([foldernames[0],foldernames[1],foldernames[2],foldernames[3],foldernames[4]])
bx.legend([foldernames[0],foldernames[1],foldernames[2],foldernames[3],foldernames[4]])

Ut.plotparams(ax,'y','time',[550,600],[-2.0,4.0],False,True,1.0,True,True,5.0)
Ut.plotparams(bx,'y','time',[550,600],[-1.25,2.75],True,True,1.0,True,True,5.0)

fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\Domain.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)
fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\Domain.png', dpi=150,format='png',bbox_inches='tight', pad_inches=0.1)





#def plotparams(plotloc,ylabel,xlabel,xlim,ylim,isx = False, yformatter = False, yformval = 0.0, isy= True, xformatter = False, xformval = 0.0):    

