import os
import numpy as np
import matplotlib.pyplot as plt
import Atul_myplt as myplt
import Utils as Ut
import math
import Atul_read_vicar_data as rvd
from matplotlib.patches import Ellipse,Arc
# import svgutils.compose as sc
# from IPython.display import SVG




#65100,65250,65300
#O 4-62300, 4.5-94000, 5-91050, 6-78150, 6.5-66400 7-85800, 8- 84700
#D 4-61050, 4.5-105500, 5-70400, 6-71300, 6.5-68900, 7-85600
pi = 4.0*math.atan(1.0)
Nbody = 2

Ur = [4,4.5,5]      #3.5
Ofilename = [62300, 94000, 91050]   #137850
Dfilename = [61050, 105500, 85200]  #129200
lcolor = ['xkcd:black','xkcd:purple','xkcd:green']  #'xkcd:gr
lstyle = ['solid','dashed','dotted']
legend = ['$U_R=4$','$U_R=4.5$','$U_R=5$']
plts1 = [[],[]]
D = 0
O = 1

figW, figH = myplt.figsize_mm(figSWmm=100, figAR=0.8,
                              nRows=3, nCols=2)

fig, axs = plt.subplots(figsize=(figW, figH),
                        nrows=3, ncols=2,
                        sharex=False,squeeze=False)

myplt.margins_adjust(fig, figW, figH, left=10, right=10, top=10, bottom=10,
                   wspace=4,hspace=2)




folder_path = 'I:\\NIHAR\\Tandem_D_O\\D_O\\'
for pl in range(len(Ur)):
    marker_path = folder_path + 'D_O_' + str(Ur[pl]) + '\\marker\\'
    plt_path = folder_path + 'D_O_' + str(Ur[pl]) + '\\plt\\'
    body = 0
    while body < Nbody:
        #markerfile name and path
        print(pl,"{:07d}".format(Ofilename[pl] if body==O else Dfilename[pl]))
        markerfilename = 'marker.' + "{:07d}".format(Ofilename[pl] if body==O else Dfilename[pl]) +  '.dat'
        pltfilename = 'q' + "{:07d}".format(Ofilename[pl] if body==O else Dfilename[pl]) + '.plt'
    
        pltfile = os.path.join(plt_path,pltfilename)
        markerfile = os.path.join(marker_path,markerfilename)
    
        xc,yc,qqf = rvd.read_binary(pltfile)
        p0 = qqf[0,0,0,2]
        print(p0)
    
        #Identify the line numbers for each body start and number of points to be read for each body
        Npoint = []
        lineidx = []
        k = 1
        with open(markerfile,'r') as f:
            while(1):
                line = f.readline()
                if line == '':
                    break
                line = line.split()
                if(line[0] == 'ZONE'):
                    Npoint.append(int(int(line[2])/3))
                    lineidx.append(k)
            
                k = k+1
        
        
        print(k)
        #read the file data for each body
    
        with open(markerfile,'r') as f:
            data = np.loadtxt(f,dtype=float,skiprows=lineidx[body],max_rows=int(Npoint[body]))
            x = data[:,0]
            y = data[:,1]
            p = data[:,-1]
        
        yc = 0.5*(max(y)+min(y))
        if body == O:
            xc = 0.5*(max(x)+min(x)) 
            r1 = (max(x)-min(x))/2
            r2 = (max(y)-min(y))/2
            #print(xc,yc,r1,r2,pi)    
            r = 0.5*(r1+r2) 
            
        if body == D:
           xc = min(x)
           r = max(x)-min(x)
           
        C = (x-xc)/r
        S = (y-yc)/r
        del x,y,xc,yc,r
         
        
        theta = np.zeros([len(S),2],dtype=float)
    
        i = 0
        while i < len(theta):
            if body == D:
                if abs(C[i] - 0.0) < 1e-5 and S[i] >= 0:
                    theta[i,0] = pi/2
                elif abs(C[i]-0.0) < 1e-5 and S[i] < 0:
                    theta[i,0] = -pi/2
                else:
                    theta[i,0] = math.atan(S[i]/C[i])
                    if C[i] < 0 and S[i] < 0:
                        theta[i,0] = theta[i,0] - pi
                    elif C[i] < 0 and S[i] > 0:
                        theta[i,0] = theta[i,0] + pi
            
            if body == O:
                theta[i,0] = math.atan(S[i]/C[i])
                if C[i] < 0 and S[i] < 0:
                    theta[i,0] = theta[i,0] - pi
                elif C[i] < 0 and S[i] > 0:
                    theta[i,0] = theta[i,0] + pi
            
            theta[i,1] = 2*(p[i]-p0)
            i = i+1
            
        theta = theta[theta[:,0].argsort()]
       
        
        ax = axs[0,body]
        plts, = ax.plot(theta[:,0]*180/pi,theta[:,1], label=legend[pl],c=lcolor[pl], ls=lstyle[pl])
        plts1[body].append(plts)
            
        body = body+1


ax = axs[0,0]
ax.legend(handles=plts1[0])
c = Arc(xy=(0, -4.7), width=45, height=1.8, theta1=-90,theta2=90,
            edgecolor='k', fc='None', lw=1)
ax.add_patch(c)
c1 = Arc(xy=(0, -4.7), width=9, height=0.3,angle= 0,theta1 = 0,theta2 = 1.75,
            edgecolor='k', fc='None', lw=0.3)
ax.add_patch(c1)
ax.plot([0,0],[-5.6,-3.8],'k-',lw=1)
ax.plot([0,23],[-4.7,-4.7],'k-',lw=0.5)
ax.plot([0,19],[-4.7,-4.7+0.7],'k-',lw=0.5)
ax.text(0.57, 0.22, '$\u03B8$', verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes, color='black', fontsize=8)
Ut.plotparams(ax,'$C_p$','$\u03B8$',[-90,90],[-6,1.0],False,True,1.0,True,True,45)
ax.text(-0.08, 0.93, '(a)',
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=12)

ax = axs[0,1]
ax.legend(handles=plts1[1])
c = Ellipse(xy=(0, 0.2), width=90, height=0.5, 
            edgecolor='k', fc='None', lw=1)
ax.add_patch(c)
c1 = Arc(xy=(0, 0.2), width=15, height=0.1,angle= 0,theta1 = 0,theta2 = 0.25,
            edgecolor='k', fc='None', lw=0.3)
ax.add_patch(c1)
ax.plot([0,45],[0.2,0.2],'k-',lw=0.5)
ax.plot([0,0],[0.2,0.2+0.25],'k-',lw=0.5)
ax.plot([0,39],[0.2,0.2+0.2],'k-',lw=0.5)
ax.text(0.56, 0.89, '$\u03B8$', verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes, color='black', fontsize=8)
Ut.plotparams(ax,'$C_p$','$\u03B8$',[-180,180],[-1.5,0.5],False,True,0.4,False,True,45)
ax.text(-0.08, 0.93, '(b)',
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=12)


Ur = [6,6.5,7]
Ofilename = [111900, 67050, 116400]
Dfilename = [120650, 68900, 116350]
lcolor = ['xkcd:black','xkcd:purple','xkcd:green']
lstyle = ['solid','dashed','dotted']
legend = ['$U_R=6$','$U_R=6.5$','$U_R=7$']
plts1 = [[],[]]
D = 0
O = 1

folder_path = 'I:\\NIHAR\\Tandem_D_O\\D_O\\'
for pl in range(len(Ur)):
    marker_path = folder_path + 'D_O_' + str(Ur[pl]) + '\\marker\\'
    plt_path = folder_path + 'D_O_' + str(Ur[pl]) + '\\plt\\'
    body = 0
    while body < Nbody:
        #markerfile name and path
        print(pl,"{:07d}".format(Ofilename[pl] if body==O else Dfilename[pl]))
        markerfilename = 'marker.' + "{:07d}".format(Ofilename[pl] if body==O else Dfilename[pl]) +  '.dat'
        pltfilename = 'q' + "{:07d}".format(Ofilename[pl] if body==O else Dfilename[pl]) + '.plt'
    
        pltfile = os.path.join(plt_path,pltfilename)
        markerfile = os.path.join(marker_path,markerfilename)
    
        xc,yc,qqf = rvd.read_binary(pltfile)
        p0 = qqf[0,0,0,2]
        print(p0)
    
        #Identify the line numbers for each body start and number of points to be read for each body
        Npoint = []
        lineidx = []
        k = 1
        with open(markerfile,'r') as f:
            while(1):
                line = f.readline()
                if line == '':
                    break
                line = line.split()
                if(line[0] == 'ZONE'):
                    #print(line)
                    Npoint.append(int(int(line[2])/3))
                    lineidx.append(k)
            
                k = k+1
                
    
        #read the file data for each body
    
        with open(markerfile,'r') as f:
            data = np.loadtxt(f,dtype=float,skiprows=lineidx[body],max_rows=int(Npoint[body]))
            x = data[:,0]
            y = data[:,1]
            p = data[:,-1]
        
        yc = 0.5*(max(y)+min(y))
        if body == O:
            xc = 0.5*(max(x)+min(x)) 
            r1 = (max(x)-min(x))/2
            r2 = (max(y)-min(y))/2
            #print(xc,yc,r1,r2,pi)    
            r = 0.5*(r1+r2) 
            
        if body == D:
           xc = min(x)
           r = max(x)-min(x)
           
        C = (x-xc)/r
        S = (y-yc)/r
        del x,y,xc,yc,r
         
        
        theta = np.zeros([len(S),2],dtype=float)
    
        i = 0
        while i < len(theta):
            if body == D:
                if abs(C[i] - 0.0) < 1e-5 and S[i] >= 0:
                    theta[i,0] = pi/2
                elif abs(C[i]-0.0) < 1e-5 and S[i] < 0:
                    theta[i,0] = -pi/2
                else:
                    theta[i,0] = math.atan(S[i]/C[i])
                    if C[i] < 0 and S[i] < 0:
                        theta[i,0] = theta[i,0] - pi
                    elif C[i] < 0 and S[i] > 0:
                        theta[i,0] = theta[i,0] + pi
            
            if body == O:
                theta[i,0] = math.atan(S[i]/C[i])
                if C[i] < 0 and S[i] < 0:
                    theta[i,0] = theta[i,0] - pi
                elif C[i] < 0 and S[i] > 0:
                    theta[i,0] = theta[i,0] + pi
            
            theta[i,1] = 2*(p[i]-p0)
            i = i+1
            
        theta = theta[theta[:,0].argsort()]
       
        
        ax = axs[1,body]
        plts, = ax.plot(theta[:,0]*180/pi,theta[:,1], label=legend[pl],c=lcolor[pl], ls=lstyle[pl])
        plts1[body].append(plts)
            
        body = body+1


ax = axs[1,0]
ax.legend(handles=plts1[0])
c = Arc(xy=(0, -3.0), width=50, height=1.5, theta1=-90,theta2=90,
            edgecolor='k', fc='None', lw=1)
ax.add_patch(c)
c1 = Arc(xy=(0, -3.0), width=9, height=0.3,angle= 0,theta1 = 0,theta2 = 1.75,
            edgecolor='k', fc='None', lw=0.3)
ax.add_patch(c1)
ax.plot([0,0],[-3.75,-2.25],'k-',lw=1)
ax.plot([0,25],[-3.0,-3.0],'k-',lw=0.5)
ax.plot([0,25],[-3.0,-3.0+0.5],'k-',lw=0.5)
ax.text(0.58, 0.25, '$\u03B8$', verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes, color='black', fontsize=8)
Ut.plotparams(ax,'$C_p$','$\u03B8$',[-90,90],[-4,1.0],False,True,1.0,True,True,45)
ax.text(-0.08, 0.93, '(c)',
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=12)

ax = axs[1,1]
ax.legend(handles=plts1[1])
c = Ellipse(xy=(0, 0.0), width=90, height=0.72, 
            edgecolor='k', fc='None', lw=1)
ax.add_patch(c)
c1 = Arc(xy=(0, 0.0), width=15, height=0.1,angle= 0,theta1 = 0,theta2 = 0.25,
            edgecolor='k', fc='None', lw=0.3)
ax.add_patch(c1)
ax.plot([0,45],[0.0,0.0],'k-',lw=0.5)
ax.plot([0,0],[0.0,0.0+0.36],'k-',lw=0.5)
ax.plot([0,40],[0.0,0.0+0.2],'k-',lw=0.5)
ax.text(0.58, 0.87, '$\u03B8$', verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes, color='black', fontsize=8)

Ut.plotparams(ax,'$C_p$','$\u03B8$',[-180,180],[-2.5,0.5],False,True,0.4,False,True,45)
ax.text(-0.08, 0.93, '(d)',
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=12)

Ur = [8.5,10,12.5]
Ofilename = [100450,100050, 80950]
Dfilename = [100500,110350, 80250]
lcolor = ['xkcd:green','xkcd:black','xkcd:purple']
lstyle = ['dotted','solid','dashed']
legend = ['$U_R=8.5$','$U_R=10$','$U_R=12.5$']
plts1 = [[],[]]
D = 0
O = 1

folder_path = 'I:\\NIHAR\\Tandem_D_O\\D_O\\'
for pl in range(len(Ur)):
    marker_path = folder_path + 'D_O_' + str(Ur[pl]) + '\\marker\\'
    plt_path = folder_path + 'D_O_' + str(Ur[pl]) + '\\plt\\'
    body = 0
    while body < Nbody:
        #markerfile name and path
        print(pl,"{:07d}".format(Ofilename[pl] if body==O else Dfilename[pl]))
        markerfilename = 'marker.' + "{:07d}".format(Ofilename[pl] if body==O else Dfilename[pl]) +  '.dat'
        pltfilename = 'q' + "{:07d}".format(Ofilename[pl] if body==O else Dfilename[pl]) + '.plt'
    
        pltfile = os.path.join(plt_path,pltfilename)
        markerfile = os.path.join(marker_path,markerfilename)
    
        xc,yc,qqf = rvd.read_binary(pltfile)
        p0 = qqf[0,0,0,2]
        print(p0)
    
        #Identify the line numbers for each body start and number of points to be read for each body
        Npoint = []
        lineidx = []
        k = 1
        with open(markerfile,'r') as f:
            while(1):
                line = f.readline()
                if line == '':
                    break
                line = line.split()
                if(line[0] == 'ZONE'):
                    #print(line)
                    Npoint.append(int(int(line[2])/3))
                    lineidx.append(k)
            
                k = k+1
                
    
        #read the file data for each body
    
        with open(markerfile,'r') as f:
            data = np.loadtxt(f,dtype=float,skiprows=lineidx[body],max_rows=int(Npoint[body]))
            x = data[:,0]
            y = data[:,1]
            p = data[:,-1]
        
        yc = 0.5*(max(y)+min(y))
        if body == O:
            xc = 0.5*(max(x)+min(x)) 
            r1 = (max(x)-min(x))/2
            r2 = (max(y)-min(y))/2
            #print(xc,yc,r1,r2,pi)    
            r = 0.5*(r1+r2) 
            
        if body == D:
           xc = min(x)
           r = max(x)-min(x)
           
        C = (x-xc)/r
        S = (y-yc)/r
        del x,y,xc,yc,r
         
        
        theta = np.zeros([len(S),2],dtype=float)
    
        i = 0
        while i < len(theta):
            if body == D:
                if abs(C[i] - 0.0) < 1e-5 and S[i] >= 0:
                    theta[i,0] = pi/2
                elif abs(C[i]-0.0) < 1e-5 and S[i] < 0:
                    theta[i,0] = -pi/2
                else:
                    theta[i,0] = math.atan(S[i]/C[i])
                    if C[i] < 0 and S[i] < 0:
                        theta[i,0] = theta[i,0] - pi
                    elif C[i] < 0 and S[i] > 0:
                        theta[i,0] = theta[i,0] + pi
            
            if body == O:
                theta[i,0] = math.atan(S[i]/C[i])
                if C[i] < 0 and S[i] < 0:
                    theta[i,0] = theta[i,0] - pi
                elif C[i] < 0 and S[i] > 0:
                    theta[i,0] = theta[i,0] + pi
            
            theta[i,1] = 2*(p[i]-p0)
            i = i+1
            
        theta = theta[theta[:,0].argsort()]
       
        
        ax = axs[2,body]
        plts, = ax.plot(theta[:,0]*180/pi,theta[:,1], label=legend[pl],c=lcolor[pl], ls=lstyle[pl])
        plts1[body].append(plts)
            
        body = body+1


ax = axs[2,0]
ax.legend(handles=plts1[0])
c = Arc(xy=(0, -3.0), width=50, height=1.5, theta1=-90,theta2=90,
            edgecolor='k', fc='None', lw=1)
ax.add_patch(c)
c1 = Arc(xy=(0, -3.0), width=9, height=0.3,angle= 0,theta1 = 0,theta2 = 1.75,
            edgecolor='k', fc='None', lw=0.3)
ax.add_patch(c1)
ax.plot([0,0],[-3.75,-2.25],'k-',lw=1)
ax.plot([0,25],[-3.0,-3.0],'k-',lw=0.5)
ax.plot([0,25],[-3.0,-3.0+0.5],'k-',lw=0.5)
ax.text(0.58, 0.25, '$\u03B8$', verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes, color='black', fontsize=8)
Ut.plotparams(ax,'$C_p$','$\u03B8$ (degrees)',[-90,90],[-4,1.0],True,True,1.0,True,True,45)
ax.text(-0.08, 0.93, '(e)',
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=12)

ax = axs[2,1]
ax.legend(handles=plts1[1])
c = Ellipse(xy=(0, 0.75), width=90, height=0.95, 
            edgecolor='k', fc='None', lw=1)
ax.add_patch(c)
c1 = Arc(xy=(0, 0.75), width=15, height=0.1,angle= 0,theta1 = 0,theta2 = 0.25,
            edgecolor='k', fc='None', lw=0.3)
ax.add_patch(c1)
ax.plot([0,45],[0.75,0.75],'k-',lw=0.5)
ax.plot([0,0],[0.75,0.75+0.475],'k-',lw=0.5)
ax.plot([0,40],[0.75,0.75+0.33],'k-',lw=0.5)
ax.text(0.58, 0.82, '$\u03B8$', verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes, color='black', fontsize=8)

Ut.plotparams(ax,'$C_p$','$\u03B8$ (degrees)',[-180,180],[-2.0,1.5],True,True,0.4,False,True,45)
ax.text(-0.08, 0.93, '(f)',
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=12)


fig.savefig('E:\\BACKUP_NIHAR_1_OCT_2020\\POF1\Cpt.pdf', dpi=150,format='pdf',bbox_inches='tight', pad_inches=0.1)
        
    