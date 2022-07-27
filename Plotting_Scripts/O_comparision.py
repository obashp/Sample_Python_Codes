import os
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import gridspec
import matplotlib.lines as mlines
import Atul_myplt as myplt
import Utils as Ut

majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1)


 


figW, figH = myplt.figsize_mm(figSWmm=200, figAR=0.4,
                              nRows=2, nCols=1)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=2, ncols=1,
                        sharex=True,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=10,
                  wspace=2,hspace=4)


    
folder_path0 = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\'

folder_path1 = 'D_O\\'
Ampfile = os.path.join(folder_path0+folder_path1,'2y_AC_Ur.dat')
Ur0 = Ut.GetUrspan(folder_path0+folder_path1)
# print(Ur0)
AmaxDO = np.zeros(len(Ur0),dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    AmaxDO[:] = data[0:len(Ur0),2]
    
folder_path1 = 'Cylinder_2DOF_VIV\\M_2.5pi'
Ampfile = os.path.join(folder_path0+folder_path1,'1y_AC_Ur.dat')
Ur3 = Ut.GetUrspan(folder_path0+folder_path1)
AmaxOB = np.zeros(len(Ur3),dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    AmaxOB[:] = data[0:len(Ur3),2]


folder_path1 = 'O_O\\'
Ampfile = os.path.join(folder_path0+folder_path1,'O_downstream.dat')   
Ur1 = np.zeros(52,dtype=float)
AmaxOOM = np.zeros(52,dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur1 = data[:,0]
    AmaxOOM[:] = data[:,1]

folder_path1 = 'Cylinder_2DOF_VIV\\Bao_Validation'
Ampfile = os.path.join(folder_path0+folder_path1,'Ym_Ur_2.dat')
Ur2 = np.zeros(10,dtype=float)
AmaxOOB = np.zeros(10,dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur2 = data[:,0]
    AmaxOOB[:] = data[:,1]
    
folder_path1 = 'Cylinder_2DOF_VIV'
Ampfile = os.path.join(folder_path0+folder_path1,'Cylinder.dat')
with open(Ampfile,'r') as f:
    data = np.loadtxt(f)
    Ur = data[:,0]
    AmaxOOW = data[:,1]

folder_path1 = 'Cylinder_2DOF_VIV'
Ampfile = os.path.join(folder_path0+folder_path1,'Jester.dat')
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,delimiter=',')
    Ur5 = data[:,0]
    AmaxOOW1 = data[:,1]

    



plotformat = ['bv-.','r+--','kv-','c<--' ]

xlim = [Ur0[0],Ur0[len(Ur0)-1]]
ax = axs[1,0]
ax.plot(Ur0,AmaxDO,ls='solid',marker='v',color='xkcd:black')
ax.plot(Ur3,AmaxOB,ls='dashed',marker='+',color='xkcd:red')
#ax.plot(Ur,AmaxOOW,ls='none',marker='o',color='xkcd:blue',mfc='none')
ax.plot(Ur1,AmaxOOM,ls='dashdot',marker='+',color='xkcd:purple',mfc='none')
ax.plot(Ur2,AmaxOOB,ls='dashed',marker='o',color='xkcd:green',mfc='none')
ax.plot(Ur5,AmaxOOW1,ls='dotted',marker='H',color='xkcd:orange',mfc='none')
ax.legend(['$D - O$ ($M = 7.85$)','Isolated $O$ ($M = 7.85$)','$O - O$ ($M = 10, Re = 100$), Prasanth et. al','$O - O$ ($M=2, Re=150$) Bao et. al','$O-O$ ($Re = 1000$), Jester et. al'],fontsize = 7.5,loc='upper left',borderpad=0.0,edgecolor='1.0')
Ut.plotparams(ax,'$A_{max}$','$U_R$',xlim,[0.0,1.55],True,True,0.3,True,True,1.0)
# ax.plot(Ur,phasepi,linestyle='--',color='grey')  
# ax.plot(Ur,phasepi*0.5,linestyle='--',color='grey')  
# ax.plot(Ur,phasepi*0,linestyle='--',color='grey')  


folder_path1 = 'D_O\\'
Ampfile = os.path.join(folder_path0+folder_path1,'1y_AC_Ur.dat')
Ur0 = Ut.GetUrspan(folder_path0+folder_path1)
AmaxDO = np.zeros(len(Ur0),dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    AmaxDO[:] = data[0:len(Ur0),2]


folder_path1 = 'IsolatedD\\'
Ampfile = os.path.join(folder_path0+folder_path1,'1y_A_Sty_Ur.dat')  
Ur1 = Ut.GetUrspan(folder_path0+folder_path1)
AmaxD = np.zeros(len(Ur1),dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur1 = data[:,0]
    AmaxD[:] = data[:,1]
    
Ampfile = os.path.join(folder_path0+folder_path1,'A_D_Thompson.dat')   
Ur2 = np.zeros(77,dtype=float)
AmaxDT = np.zeros(77,dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur2 = data[:,0]
    AmaxDT[:] = data[:,1]

labels = ['(a)','(b)']
for t in range(2):
    ax = axs[t,0]
    ax.text(-0.05, 0.93, labels[t],
        verticalalignment='top', horizontalalignment='right',   
        transform=ax.transAxes,
        color='black', fontsize=15)

    
plotformat = ['bv-.','r+--','kv-','c+--' ]


ax = axs[0,0]
ax.plot(Ur0,AmaxDO,ls='solid',marker='v',color='xkcd:black')
ax.plot(Ur1,AmaxD,ls='dashed',marker='+',color='xkcd:red')
ax.plot(Ur2,AmaxDT,ls='dashed',marker='+',color='xkcd:green',mfc='none')
ax.legend(['$D - O$ (Present study)','Isolated $D$ ($Re = 100$)','Isolated $D$ ($M = 6.0$, $Re \in [1080,9000]$), Zhao et. al'],fontsize=7.5,edgecolor='1.0')
#ax.legend(['$DO - O-A_{10}$ (Present study)','$OO - O-A_{10}$ (Prasanth et.al)','$OO-O-A_{10}$ (Bao et.al)','$O-A_{10}$ (Bao et.al)'])
Ut.plotparams(ax,'$A_{max}$','$U_R$',xlim,[0.0,3.5],False,True,0.4,True,True,1.0)




     
fig.tight_layout()
fig.savefig('E:\BACKUP_NIHAR_1_OCT_2020\\POF1\\AmpcompO.pdf', dpi=100,format='pdf',bbox_inches='tight', pad_inches=0.1)


figW, figH = myplt.figsize_mm(figSWmm=150, figAR=0.4,
                              nRows=2, nCols=2)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=2, ncols=2,
                        sharex=True,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=10, right=10, top=5, bottom=5,
                  wspace=4,hspace=6)


folder_path0 = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\'

folder_path1 = 'Cylinder_2DOF_VIV\\Bao_Validation\\'


Ampfile = os.path.join(folder_path0+folder_path1,'1x_AC_Ur.dat')
Ur0 = Ut.GetUrspan(folder_path0+folder_path1)
AmaxDO = np.zeros(len(Ur0),dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    AmaxDO[:] = data[0:len(Ur0),1]

folder_path1 = 'Cylinder_2DOF_VIV\\Bao_Validation'
Ampfile = os.path.join(folder_path0+folder_path1,'Xm_Ur_1.dat')
Ur2 = np.zeros(10,dtype=float)
AmaxOOB = np.zeros(10,dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur2 = data[:,0]
    AmaxOOB[:] = data[:,1]
    
ax = axs[0,0]
ax.plot(Ur0,AmaxDO,ls='solid',marker='v',color='xkcd:black')
ax.plot(Ur2,AmaxOOB,ls='dashed',marker='+',color='xkcd:red')
ax.legend(['$O - O$ (Present study)','$O - O$ (Bao et. al)'],fontsize=8,edgecolor='1.0')
ax.set_title('Upstream cylinder - In-Line response')
Ut.plotparams(ax,'$A_{X_{max}}$','$U_R$',[2.5,12.5],[0.0,0.1],False,True,0.03,True,True,1.0)
ax.text(-0.08, 0.93, '(a)',
        verticalalignment='top', horizontalalignment='right',   
        transform=ax.transAxes,
        color='black', fontsize=15)

Ampfile = os.path.join(folder_path0+folder_path1,'2x_AC_Ur.dat')
Ur0 = Ut.GetUrspan(folder_path0+folder_path1)
AmaxDO = np.zeros(len(Ur0),dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    AmaxDO[:] = data[0:len(Ur0),1]

folder_path1 = 'Cylinder_2DOF_VIV\\Bao_Validation'
Ampfile = os.path.join(folder_path0+folder_path1,'Xm_Ur_2.dat')
Ur2 = np.zeros(10,dtype=float)
AmaxOOB = np.zeros(10,dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur2 = data[:,0]
    AmaxOOB[:] = data[:,1]
    
ax = axs[0,1]
ax.plot(Ur0,AmaxDO,ls='solid',marker='v',color='xkcd:black')
ax.plot(Ur2,AmaxOOB,ls='dashed',marker='+',color='xkcd:red')
ax.legend(['$O - O$ (Present study)','$O - O$ (Bao et. al)'],fontsize=8,edgecolor='1.0')
ax.set_title('Downstream cylinder - In-line response')
Ut.plotparams(ax,'$A_{X_{max}}$','$U_R$',[2.5,12.5],[0.0,0.35],False,True,0.08,True,True,1.0)
ax.text(-0.08, 0.93, '(b)',
        verticalalignment='top', horizontalalignment='right',   
        transform=ax.transAxes,
        color='black', fontsize=15)

Ampfile = os.path.join(folder_path0+folder_path1,'1y_AC_Ur.dat')
Ur0 = Ut.GetUrspan(folder_path0+folder_path1)
AmaxDO = np.zeros(len(Ur0),dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    AmaxDO[:] = data[0:len(Ur0),1]

Ampfile = os.path.join(folder_path0+folder_path1,'Ym_Ur_1.dat')
Ur2 = np.zeros(10,dtype=float)
AmaxOOB = np.zeros(10,dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur2 = data[:,0]
    AmaxOOB[:] = data[:,1]
    
ax = axs[1,0]
ax.plot(Ur0,AmaxDO,ls='solid',marker='v',color='xkcd:black')
ax.plot(Ur2,AmaxOOB,ls='dashed',marker='+',color='xkcd:red')
# ax.legend(['$O - O$ (Present study)','$O - O$ (Bao et. al)'],fontsize=8,edgecolor='1.0')
ax.set_title('Upstream cylinder - Transverse response')
Ut.plotparams(ax,'$A_{Y_{max}}$','$U_R$',[2.5,12.5],[0.0,0.75],True,True,0.15,True,True,1.0)
ax.text(-0.08, 0.93, '(c)',
        verticalalignment='top', horizontalalignment='right',   
        transform=ax.transAxes,
        color='black', fontsize=15)

folder_path1 = 'Cylinder_2DOF_VIV\\Bao_Validation\\'
Ampfile = os.path.join(folder_path0+folder_path1,'2y_AC_Ur.dat')
Ur0 = Ut.GetUrspan(folder_path0+folder_path1)
AmaxDO = np.zeros(len(Ur0),dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    AmaxDO[:] = data[0:len(Ur0),1]

folder_path1 = 'Cylinder_2DOF_VIV\\Bao_Validation'
Ampfile = os.path.join(folder_path0+folder_path1,'Ym_Ur_2.dat')
Ur2 = np.zeros(10,dtype=float)
AmaxOOB = np.zeros(10,dtype=float)
with open(Ampfile,'r') as f:
    data = np.loadtxt(f,skiprows=1)
    Ur2 = data[:,0]
    AmaxOOB[:] = data[:,1]
    
ax = axs[1,1]
ax.plot(Ur0,AmaxDO,ls='solid',marker='v',color='xkcd:black')
ax.plot(Ur2,AmaxOOB,ls='dashed',marker='+',color='xkcd:red')
# ax.legend(['$O - O$ (Present study)','$O - O$ (Bao et. al)'],fontsize=8,edgecolor='1.0')
ax.set_title('Downstream cylinder - Transverse response')
Ut.plotparams(ax,'$A_{Y_{max}}$','$U_R$',[2.5,12.5],[0.0,1.25],True,True,0.2,True,True,1.0)
ax.text(-0.08, 0.93, '(d)',
        verticalalignment='top', horizontalalignment='right',   
        transform=ax.transAxes,
        color='black', fontsize=15)

fig.tight_layout()
fig.savefig('E:\BACKUP_NIHAR_1_OCT_2020\\POF1\\O_comparison.pdf', dpi=100,format='pdf',bbox_inches='tight', pad_inches=0.1)
    