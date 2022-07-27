import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Atul_myplt as myplt
import Utils as Ut
from matplotlib.patches import Ellipse
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)


folder_path = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\D_O\\'
folder_path1 = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\IsolatedD\\'
folder_path2 = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\Cylinder_2DOF_VIV\\M_2.5pi'

body = '1'
dire = 'y'
St = 0.1680


A = '_AC_Ur.dat'
Ampfile1 = os.path.join(folder_path,'1y'+A)
Ampfile2 = os.path.join(folder_path,'2y'+A)
Ampfile3 = os.path.join(folder_path1,'1y'+A)
Ampfile4 = os.path.join(folder_path2,'1y'+A)

Ur = Ut.GetUrspan(folder_path)
nUr = len(Ur)
print(nUr)



figW, figH = myplt.figsize_mm(figSWmm=150, figAR=0.8,
                              nRows=1, nCols=1)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=1, ncols=1,
                        sharex=True,squeeze=False)

myplt.margins_adjust(fig, figW, figH, left=10, right=10, top=5, bottom=5,
                   wspace=0,hspace=2)


with open(Ampfile1,'r') as f:
    data = np.loadtxt(Ampfile1,skiprows=1)
    Ur = data[0:nUr,0]
    Sty1 = data[0:nUr,7]

with open(Ampfile2,'r') as f:
    data = np.loadtxt(Ampfile2,skiprows=1)    
    Sty2 = data[0:nUr,7]
    
with open(Ampfile3,'r') as f:
    data = np.loadtxt(Ampfile3,skiprows=1)
    Ur1 = data[0:nUr,0]
    Sty3 = data[0:nUr,10]
    
with open(Ampfile4,'r') as f:
    data = np.loadtxt(Ampfile4,skiprows=1)
    Ur2 = data[0:nUr,0]
    Sty4 = data[0:nUr,7]

Fn = 1/Ur;    

axs[0,0].plot(Ur,Sty1,linestyle='none',color='xkcd:black',mfc='None',marker='v')
axs[0,0].plot(Ur,Sty2,linestyle='none',color='xkcd:red',mfc='None',marker='o')
#axs[0,0].plot([min(Ur),max(Ur)],[0.168,0.168])
axs[0,0].plot(Ur1,Sty3,linestyle='none',color='xkcd:blue',marker='s',mfc='None')
axs[0,0].plot(Ur2,Sty4,linestyle='none',color='xkcd:purple',marker='H',mfc='None')
axs[0,0].plot(Ur,Fn,linestyle='dotted',color='xkcd:gray')

axs[0,0].set_xlim([1,15])
axs[0,0].set_ylim([1/17,0.28])

xmajorLocator = MultipleLocator(2)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(0.5)
axs[0,0].xaxis.set_minor_locator(xminorLocator)
axs[0,0].xaxis.set_major_locator(xmajorLocator)

ymajorLocator = MultipleLocator(0.04)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator = MultipleLocator(0.01)
axs[0,0].yaxis.set_minor_locator(yminorLocator)
axs[0,0].yaxis.set_major_locator(ymajorLocator)

axs[0,0].set_ylabel('$f_y$',fontsize=12,color='black')
axs[0,0].set_xlabel('$U_R$',fontsize=12,color='black')

axs[0,0].tick_params(which='both',left=1, top=1, right=1,
                    bottom=1, labelleft=1, labeltop=0,
                    labelright=0, labelbottom=1, width=0.5, direction="in")

axs[0,0].legend(['D cylinder (D-O)','O cylinder (D-O)','D cylinder (Isolated)','O cylinder (Isolated)','$f_N = 1/U_R$'])

fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\Sty_Ur.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)



