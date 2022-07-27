import os
import numpy as np
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

def streams(ax,xx,yy,u,v,speed,base_map=False):
 x = np.linspace(xx.min(), xx.max(), 50)
 y = np.linspace(yy.min(), yy.max(), 50)

 xi, yi = np.meshgrid(x,y)

 #then, interpolate your data onto this grid:

 px = xx.flatten()
 py = yy.flatten()
 pu = u.flatten()
 pv = v.flatten()
 pspeed = speed.flatten()

 gu = griddata(zip(px,py), pu, (xi,yi))
 gv = griddata(zip(px,py), pv, (xi,yi))
 gspeed = griddata(zip(px,py), pspeed, (xi,yi))

 lw = 6*gspeed/np.nanmax(gspeed)
 #now, you can use x, y, gu, gv and gspeed in streamplot:

 if base_map:
    xx,yy = ax(xx,yy)
    xi,yi = ax(xi,yi)

 ax.contour(xx,yy,speed, colors='k', alpha=0.4)
 ax.plot(xx,yy,'-k',alpha=0.3)
 ax.plot(xx.T,yy.T,'-k',alpha=0.3)
 ax.plot(xi,yi,'-b',alpha=0.1)
 ax.plot(xi.T,yi.T,'-b',alpha=0.1)
 c = ax.streamplot(x,y,gu,gv, density=2,
    linewidth=lw, color=gspeed, cmap=plt.cm.jet)