import os
import numpy as np
import matplotlib.pyplot as plt
import Atul_myplt as myplt
import Utils as ut



folder_path0 = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\D_O'
Urfile = os.path.join(folder_path0,'Urlist.dat')
with open(Urfile, 'r') as f:
    Ur = np.loadtxt(f)
    
Ur.sort()


Urs = [3]
Uri = np.zeros(len(Urs),dtype = int)

j = 0
for i in range(len(Ur)):
    if Ur[i] == Urs[j]:
        Uri[j] = i;
        j = j+1;
    if j==len(Urs):
        break


NI = len(Ur)
NJ = 301  
NUMVAR = 2  
xc = np.zeros((NI,1),dtype=float)
yc = np.zeros((NJ,1),dtype=float)
qq = np.zeros((NJ,NI,NUMVAR),dtype=float)

#St0 = 0.128174 # Re = 50

St0 = [0.181736,0.168067]#St for D_section, W=8D at Re=100
print(St0)

St1 = np.zeros((2,NI),dtype = float)

Power = [-6,2]
StD = [0.181736, 0.181736]
StO = [0.168067, 0.168067]
fN = [1.0,1.0]

for i in range(NI):
        #St1[i] = St0*Ur[i]*np.sqrt(1.5) used for added mass effect
        St1[0,i] = St0[0]*Ur[i]
        St1[1,i] = St0[1]*Ur[i]
        

filename1 = os.path.join(folder_path0,'PSD_body_1y.plt')
filename2 = os.path.join(folder_path0,'PSD_body_2y.plt')


with open(filename1, 'r') as f: 
    data1 = np.loadtxt(f,skiprows=2)
with open(filename2,'r') as f:
    data2 = np.loadtxt(f,skiprows=2)
    

x1 = data1[:,0]
y1 = data1[:,1]

x2 = data2[:,0]
y2 = data2[:,1]

yll = y1.min();  yul = y1.max()

xc = Ur 
yc = np.linspace(yll, yul, NJ) 

    
X,Y = np.meshgrid(xc,yc)

count = 0
for i in range(NI):
    for j in range(NJ):
                qq[j,i,0] = data1[count,2]
                qq[j,i,1] = data2[count,2]
                count = count + 1


qq[:,:,0] = np.minimum(qq[:,:,0],2)
qq[:,:,0] = np.maximum(qq[:,:,0],0)


fq0 = qq[:,:,0]
fq0 = ut.safe_ln(fq0)

fq1 = qq[:,:,1]
fq1 = ut.safe_ln(fq1)


rowNum = 1#int(len(Uri)/2)  # num of rows of the subplots
colNum = 1 # num of colums of the subplots

figW, figH = myplt.figsize_mm(figSWmm=260, figAR=0.4,
                              nRows=rowNum, nCols=colNum)



fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=rowNum, ncols=colNum,                                                                                                                                                                                                             
                        sharex=False,squeeze=False)  # sharing the x-label

myplt.margins_adjust(fig, figW, figH, left=20, right=20, top=5, bottom=10,
                    wspace=2,hspace=4)



k = 0
for j in range(colNum):
    for i in range(rowNum):
        ax = axs[i,j]
        ax.plot(fN,Power,color ='0.8',ls = '-')
        ax.plot([2*fN[0],2*fN[1]],Power,color ='0.8',ls = '--')
        ax.plot([3*fN[0],3*fN[1]],Power,color ='0.8',ls = '--')
        ax.plot([4*fN[0],4*fN[1]],Power,color ='0.8',ls = '--')
        ax.plot([5*fN[0],5*fN[1]],Power,color ='0.8',ls = '--')
        ax.plot([Urs[k]*StD[0],Urs[k]*StD[1]],Power,color ='0.5',ls = '-.')
        ax.plot([Urs[k]*StO[0],Urs[k]*StO[1]],Power,color ='0.5',ls = ':')
        if(Urs[k] >= 6.5):
            ax.plot([3*Urs[k]*StD[0],3*Urs[k]*StD[1]],Power,color ='0.5',ls = '-.')
            ax.plot([3*Urs[k]*StO[0],3*Urs[k]*StO[1]],Power,color ='0.5',ls = ':')
            ax.plot([0.5*Urs[k]*StD[0],0.5*Urs[k]*StD[1]],Power,color ='0.5',ls = '-.')
            ax.plot([0.5*Urs[k]*StO[0],0.5*Urs[k]*StO[1]],Power,color ='0.5',ls = ':')
            ax.plot([1.4*fN[0],1.4*fN[1]],Power,color =(0.5,0.4,0.1),ls = '-')
            ax.plot([Urs[k]*StD[0],Urs[k]*StD[1]],Power,color ='0.3',ls = '-.')
            ax.plot([Urs[k]*StO[0],Urs[k]*StO[1]],Power,color ='0.3',ls = ':')
        #if(Urs[k] > 6.5):
           # ax.plot([2.5*fN[0],2.5*fN[1]],Power,color ='0.7',ls = '--')
            # ax.plot([1.5*Urs[k]*StD[0],1.5*Urs[k]*StD[1]],Power,color ='0.5',ls = '-.')
            # ax.plot([1.5*Urs[k]*StO[0],1.5*Urs[k]*StO[1]],Power,color ='0.5',ls = ':')
            # ax.plot([0.25*Urs[k]*StD[0],0.25*Urs[k]*StD[1]],Power,color ='0.5',ls = '-.')
            # ax.plot([0.25*Urs[k]*StO[0],0.25*Urs[k]*StO[1]],Power,color ='0.5',ls = ':')
        if(Urs[k] < 6.5 and Urs[k] > 4.5):
            ax.plot([0.5*Urs[k]*StD[0],0.5*Urs[k]*StD[1]],Power,color ='0.5',ls = '-.')
            ax.plot([0.5*Urs[k]*StO[0],0.5*Urs[k]*StO[1]],Power,color ='0.5',ls = ':')
        
            # ax.plot([2.0*Urs[k]*StD[0],2.0*Urs[k]*StD[1]],Power,color ='0.5',ls = '-.')
            # ax.plot([2.0*Urs[k]*StO[0],2.0*Urs[k]*StO[1]],Power,color ='0.5',ls = ':')
        maxfD = np.argmax(fq0[:,Uri[k]])#.index(max(fq0[:,Uri[k]]))
        maxfO = np.argmax(fq1[:,Uri[k]])#.index(max(fq1[:,Uri[k]]))
        maxfD = [yc[maxfD],yc[maxfD]]
        maxfO = [yc[maxfO],yc[maxfO]]
        ax.plot([maxfD,maxfD],Power,color =(0.5,0.0,0.5),ls = '-')
        ax.plot([maxfO,maxfO],Power,color =(0.0,0.5,0.5),ls = '-.')
        ax.plot(yc,fq0[:,Uri[k]],'b',label='D-Cy')
        ax.plot(yc,fq1[:,Uri[k]],'k',label='O-Cy')
        title = '$U_R = ' + str(Urs[k]) + '$'
        ax.set_title(title)
        isx = True if (i == rowNum-1) else False
        isy = True if (j == 0) else False
        # isx = False, yformatter = False, yformval = 0.0, isy= True, xformatter = False, xformval = 0.0
        ut.plotparams(ax,'$PSD$','$f$',[0.0,6.0],[-6.0,0.5],isx, True, 1.0, isy, True,0.5)
        plotidx = '(' + chr(97+k) + ')'
        ax.text(-0.05, 0.93, plotidx,
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=12)
        k = k+1
        
lines, labels = fig.axes[-1].get_legend_handles_labels()
fig.legend(lines,labels,loc='upper right')

#print(type(ax3))
# ax1.contourf(X,Y,fq0,50,cmap=plt.cm.afmhot_r) 
# ax2.contourf(X,Y,fq1,50,cmap=plt.cm.afmhot_r) 
# ax3.contourf(X,Y,fq2,100,cmap=plt.cm.afmhot_r) 
# ax3.plot(Ur,St1[0], linestyle='--', label = '', color = 'b')
# ax3.plot(Ur,St1[1], linestyle='--', label = '', color = 'b')
# ax3.plot(Ur,fN, linestyle='--', label = '', color = 'b')

fig.tight_layout()
fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\\frq_Comp1.pdf', dpi=300,format='pdf',bbox_inches='tight', pad_inches=0.1)


    
    
