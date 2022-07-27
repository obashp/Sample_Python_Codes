import os
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import Atul_myplt as myplt
import Utils as ut

folder_path = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\D_O\\'
Ur = ut.GetUrspan(folder_path)


rowNum = 2
colNum = 3

figW, figH = myplt.figsize_mm(figSWmm=100, figAR=0.8,
                              nRows=rowNum, nCols=colNum)



fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=rowNum, ncols=colNum,                                                                                                                                                                                                             
                        sharex=False,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=8,
                    wspace=6,hspace=6)

Urs = [3,3.5,4]
 

y = ['$y$','$C_y$','$C_v$','$v_y$','$C_x$','$x$','$v_x$']

labels = [['(aD)','(bD)','(cD)'],['(aO)','(bO)','(cO)']]
print(labels[0][1])


for i in range(3):
    for j in range(2):
        ax = axs[j,i]
        xlim = [-2.2,2.2] if j == 0 else [-2.5,2.5]
        ylim = [-0.7,0.7] if j == 0 else [-0.15,0.15]
        ut.plotaxes(ax,xlim[0],xlim[1],ylim[0],ylim[1])
        isy = True
        xfval = 0.8 
        yfval = 0.2 if j == 0 else 0.1 
        isx = True
        
        ut.plotparams(ax,y[0] if i == 0 else '',y[1] if j == 1 else '',xlim,ylim,isx,True,yfval,isy,True,xfval)
        if j == 0:
            x1, y1 = ut.phasedata('1y', folder_path, 'D_O_', Urs[i], 4, 1)
            s = '$U_R = ' + str(Urs[i]) + '$'
            ax.set_title(s)
        if j == 1:
            x1,y1 = ut.phasedata('2y', folder_path, 'D_O_', Urs[i], 4, 1)
        ax.plot(x1,y1,'k')
        ax.text(-0.08, 0.93, labels[j][i],
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=10)
        
if os.path.exists("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_1.pdf"):
  os.remove("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_1.pdf")
        
fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_1.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)

Urs = [5.5,6.5,7]

labels = [['(aD)','(bD)','(cD)'],['(aO)','(bO)','(cO)']]
figW, figH = myplt.figsize_mm(figSWmm=100, figAR=0.8,
                              nRows=rowNum, nCols=colNum)



fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=rowNum, ncols=colNum,                                                                                                                                                                                                             
                        sharex=False,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=8,
                    wspace=5,hspace=6)

for i in range(3):
    for j in range(2):
        ax = axs[j,i]
        xlim = [-0.8,0.8] if j == 0 else [-1.1,1.1]
        ylim = [-1.5,1.5] if j == 0 else [-1.3,1.3]
        ut.plotaxes(ax,xlim[0],xlim[1],ylim[0],ylim[1])
        isy = True
        xfval = 0.4 
        yfval = 0.3 if j == 0 else 0.3 
        isx = True
        
        ut.plotparams(ax,y[0] if i == 0 else '',y[1] if j == 1 else '',xlim,ylim,isx,True,yfval,isy,True,xfval)
        if j == 0:
            x1, y1 = ut.phasedata('1y', folder_path, 'D_O_', Urs[i], 4, 1)
            s = '$U_R = ' + str(Urs[i]) + '$'
            ax.set_title(s)
        if j == 1:
            x1,y1 = ut.phasedata('2y', folder_path, 'D_O_', Urs[i], 4, 1)
        ax.plot(x1,y1,'k')
        ax.text(-0.08, 0.93, labels[j][i],
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=10)
        
if os.path.exists("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_2.pdf"):
  os.remove("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_2.pdf")
        
fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_2.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)


Urs = [8,9.5,10,12.5]
rowNum = 2
colNum = 4
labels = [['(aD)','(bD)','(cD)','(dD)'],['(aO)','(bO)','(cO)','(dO)']]

figW, figH = myplt.figsize_mm(figSWmm=100, figAR=0.8,
                              nRows=rowNum, nCols=colNum)



fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=rowNum, ncols=colNum,                                                                                                                                                                                                             
                        sharex=False,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=8,
                    wspace=5,hspace=6)


for i in range(4):
    for j in range(2):
        ax = axs[j,i]
        xlim = [-1.0,1.0] if j == 0 else [-1.6,1.6]
        ylim = [-2.4,2.4] if j == 0 else [-1.2,1.2]
        ut.plotaxes(ax,xlim[0],xlim[1],ylim[0],ylim[1])
        isy = True
        xfval = 0.4 
        yfval = 0.5 if j == 0 else 0.3 
        isx = True
        
        ut.plotparams(ax,y[0] if i == 0 else '',y[1] if j == 1 else '',xlim,ylim,isx,True,yfval,isy,True,xfval)
        if j == 0:
            x1, y1 = ut.phasedata('1y', folder_path, 'D_O_', Urs[i], 4, 1)
            s = '$U_R = ' + str(Urs[i]) + '$'
            ax.set_title(s)
        if j == 1:
            x1,y1 = ut.phasedata('2y', folder_path, 'D_O_', Urs[i], 4, 1)
        ax.plot(x1,y1,'k')
        ax.text(-0.08, 0.93, labels[j][i],
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=10)
        
if os.path.exists("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_3.pdf"):
  os.remove("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_3.pdf")
        
fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_3.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)



# if os.path.exists("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_t.pdf"):
#   os.remove("E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_t.pdf")
        
# fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Cy_y_t.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)



