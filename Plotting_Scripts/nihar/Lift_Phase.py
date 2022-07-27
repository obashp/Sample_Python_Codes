import os
import numpy as np
import matplotlib.pyplot as plt
import Atul_myplt as myplt
import Utils as Ut


folder_path = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\D_O'
 

filepre = ['1y','2x','2y']
label = ['(a)','(b)','(c)','(d)','(e)']
maxC = [0.0,0.0,0.0]


A = '_AC_Ur.dat'
Ur = Ut.GetUrspan(folder_path)
nUr = len(Ur)

print(nUr)

figW, figH = myplt.figsize_mm(figSWmm=250,  figAR=0.3,
                              nRows=3, nCols=1)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=3, ncols=1,
                        sharex=True,squeeze=False)

myplt.margins_adjust(fig, figW, figH, left=10, right=10, top=5, bottom=5,
                   wspace=0,hspace=2)

Ampfile = os.path.join(folder_path,filepre[0] + A)
Ut.Cfplot(axs[0,0],Ampfile,nUr,'kv-',2)
Ut.Cfplot(axs[1,0],Ampfile,nUr,'kv-',4)
Ut.Cfplot(axs[1,0],Ampfile,nUr,'co-.',5)
Ut.Cfplot(axs[2,0],Ampfile,nUr,'kv-',6)

Ampfile = os.path.join(folder_path,filepre[2] + A)
Ut.Cfplot(axs[0,0],Ampfile,nUr,'bv-',2)
Ut.Cfplot(axs[1,0],Ampfile,nUr,'bv-',4)
Ut.Cfplot(axs[1,0],Ampfile,nUr,'rv-.',5)
Ut.Cfplot(axs[2,0],Ampfile,nUr,'bv-',6)

xlim = [min(Ur),max(Ur)]
Ut.plotparams(axs[0,0],'$A_{10}$','$U_R$',xlim,[0.0,4.0],False,True,0.5,True,True,1)
axs[0,0].legend(['$D - A_{10}$','$O - A_{10}$'])
Ut.plotparams(axs[1,0],'$C_{y-R.M.S.}$','$U_R$',xlim,[0.0,1.8],False,True,0.2,True,True,1)
axs[1,0].legend(['$D-C_{y_{R.M.S.}}}$','$D-C_{yv_{R.M.S.}}$','$O-C_{y_{R.M.S.}}$','$O-C_{yv_{R.M.S.}}$'])
Ut.plotparams(axs[2,0],'$C_{x-R.M.S.}$','$U_R$',xlim,[0.0,5.0],True,True,0.5,True,True,1)
axs[2,0].legend(['$D-C_{x_{R.M.S.}}$','$O-C_{x_{R.M.S.}}$'])



# folder_path1 = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\Cylinder_2DOF_VIV\\f1'
# Ampfile = os.path.join(folder_path1,filepre[0] + A)
# Ut.Cfplot(axs[1,0],Ampfile,nUr,'r--.',6)
# Ampfile = os.path.join(folder_path,filepre[2] + A)
# Ut.Cfplot(axs[1,0],Ampfile,nUr,'kv-',6)
# Ut.Cfplot(axs[2,0],Ampfile,nUr,'bv-',5,True)
# Ut.Cfplot(axs[1,0],Ampfile,nUr,'co-.',8)
# Ut.Cfplot(axs[2,0],Ampfile,nUr,'rv--',7,True)



for i in range(3):
    axs[i,0].text(-0.05, 0.93, label[i],
        verticalalignment='top', horizontalalignment='right',
        transform=axs[i,0].transAxes,
        color='black', fontsize=12)


fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\CyPhase.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)