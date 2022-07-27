import os
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import Atul_myplt as myplt
import Utils as ut

folder_path = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\D_O'
Ur = ut.GetUrspan(folder_path)

print(Ur)

rowNum = 3
colNum = 3

figW, figH = myplt.figsize_mm(figSWmm=150, figAR=1.0,
                              nRows=rowNum, nCols=colNum)



fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=rowNum, ncols=colNum,                                                                                                                                                                                                             
                        sharex=False,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=10,
                    wspace=2,hspace=4)

Urs = [3,3.5,4,6,6.5,7,12,12.5,13]
 


for i in range(rowNum):
    for j in range(colNum):
        ax = axs[i,j]
        k = i*rowNum+j
        print(k,Urs[k])
        ut.plotaxes(axs[i,j],-3.0,3.0,-1.0,1.0)
        if j == 0:
            if i == 2:
                ut.plotparams(axs[i,j],'y','$C_y$',[-3.0,3.0],[-1.0,1.0],True,True,0.5)
            else:
                ut.plotparams(axs[i,j],'y','$C_y$',[-3.0,3.0],[-1.0,1.0],False,True,0.5)
        else:
            if i == 2:
                ut.plotparams(axs[i,j],'y','$C_y$',[-3.0,3.0],[-1.0,1.0],True,True,0.5,False)
            else:
                ut.plotparams(axs[i,j],'y','$C_y$',[-3.0,3.0],[-1.0,1.0],False,True,0.5,False)
                
        filepath = '\\D_O_' + str(Urs[k])
        filename = '1y_phasep_body.dat'
        file = folder_path + filepath + filename
        with open(file,'r') as f:
            data = np.loadtxt(f,skiprows = 2)
        print(data)
        y = data[:,0]
        vy = data[:,1]
        ay = data[:,2]
        Cy = data[:,3]
        ax.plot(Cy,y,'k')





