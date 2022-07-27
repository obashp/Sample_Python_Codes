import os
import numpy as np
import matplotlib.pyplot as plt
import Atul_myplt as myplt
import Utils as ut

folder_path0 = 'E:\BACKUP_NIHAR_1_OCT_2020\\POF1\\'
f = ['xgrid.dat','ygrid.dat','UnstrucSurfBody1.dat','UnstrucSurfBody2.dat']

xfile = os.path.join(folder_path0,f[0])
yfile = os.path.join(folder_path0,f[1])

fig, axs = plt.subplots(figsize=(7, 4),
                        nrows=1, ncols=1,
                        sharex=False,squeeze=False)

ut.plot2Dgrid(axs[0,0],xfile,yfile,'xkcd:black solid')

B1file = os.path.join(folder_path0,'UnstrucSurfBody1.dat')
style = 'xkcd:red '+ 'solid'
ut.plot2DUnstrucSurf(axs[0,0],B1file,style) 
        
B2file = os.path.join(folder_path0,'UnstrucSurfBody2.dat')
ut.plot2DUnstrucSurf(axs[0,0],B2file,style)

fig.savefig('E:\BACKUP_NIHAR_1_OCT_2020\\POF1\\Grid_complete.pdf', dpi=300,format='pdf',bbox_inches='tight', pad_inches=0.2)

figW, figH = myplt.figsize_mm(figSWmm=140, figAR=0.9,
                              nRows=1, nCols=3)

fig1, axs1 = plt.subplots(figsize=(figW, figH),
                        nrows=1, ncols=3,
                        sharex=False,squeeze=False)

ax = axs1[0,0]
ut.plot2Dgridzoomed(ax, xfile, yfile, 'xkcd:black solid', 33, 112, 47, 146)
ut.plot2DUnstrucSurf(ax,B1file,style)
ut.plotparams(ax,'Y','X',[4.713,6.293],[9.21,10.79],True,True,0.2,True,True,0.2)
axs1[0,0].text(-0.02, 0.85, '(a)',
        verticalalignment='top', horizontalalignment='right',
        transform=axs[0,0].transAxes,
        color='black', fontsize=12)


ax = axs1[0,1]
ut.plot2Dgridzoomed(ax, xfile, yfile, 'xkcd:black solid', 261, 340, 47, 146)
ut.plot2DUnstrucSurf(ax,B2file,style)
ut.plotparams(ax,'Y','X',[9.713,11.293],[9.21,10.79],True,False,0.2,False,True,0.2)
axs1[0,1].text(0.58, 0.85, '(b)',
        verticalalignment='top', horizontalalignment='right',
        transform=axs[0,0].transAxes,
        color='black', fontsize=12)

ax = axs1[0,2]
ut.plot2Dgridzoomed(ax, xfile, yfile, 'xkcd:black solid', 400, 479, 47, 146)
ut.plotparams(ax,'Y','X',[12.493,20.25],[9.21,10.79],True,False,0.2,False,True,1.0)
axs1[0,2].text(1.18, 0.85, '(c)',
        verticalalignment='top', horizontalalignment='right',
        transform=axs[0,0].transAxes,
        color='black', fontsize=12)
fig1.savefig('E:\BACKUP_NIHAR_1_OCT_2020\\POF1\\Grid_zoomed.pdf', dpi=300,format='pdf',bbox_inches='tight', pad_inches=0.2)
  
            
        
    
