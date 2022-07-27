# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 23:01:11 2017

@author: Jisheng
"""

import numpy as np
import matplotlib
import matplotlib.pylab as plt
import matplotlib.patches as patches
import matplotlib.collections as PolyCollection
import matplotlib.transforms as mtransforms
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap
import sys
import os
import __main__
import scipy.interpolate
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numpy import NaN, Inf, arange, isscalar, asarray, array
import re
import math


###--- figsize

def figsize_mm(figSWmm=65, figAR=0.618, nRows=1, nCols=1):
    figW = nCols * figSWmm / 25.4
    figH = nRows * figSWmm * figAR / 25.4
    return figW, figH


# --- adjust margins
def margins_adjust(fig, figW, figH,
                   left=10, right=2.5, bottom=10, top=2.5,
                   wspace=15, hspace=2):
    fig.subplots_adjust(left=left / (figW * 25.4),
                        right=1 - right / (figW * 25.4),
                        bottom=bottom / (figH * 25.4),
                        top=1 - top / (figH * 25.4),
                        wspace=wspace / 25.4,
                        hspace=hspace / 25.4)


# --- save figs in different formats
def savefig(fig, figname, dpi=300):
    fig.savefig(figname + '.eps', dpi=dpi)
    fig.savefig(figname + '.pdf', dpi=dpi)
    #fig.savefig(figname + '.png', dpi=dpi)


def get_appName():
#    appName = os.path.basename(__main__.__file__).strip(".py")
    appName = os.path.basename(__main__.__file__)
    appName = os.path.splitext(appName)[0]
    return appName


# --- customise colour map
def make_colormap(colors):
    from matplotlib.colors import LinearSegmentedColormap, ColorConverter
    from numpy import sort

    z = np.array(sorted(colors.keys()))
    n = len(z)
    z1 = min(z)
    zn = max(z)
    x0 = (z - z1) / (zn - z1)

    CC = ColorConverter()
    R = []
    G = []
    B = []
    for i in range(n):
        Ci = colors[z[i]]
        if type(Ci) == str:
            RGB = CC.to_rgb(Ci)
        else:
            RGB = Ci
        R.append(RGB[0])
        G.append(RGB[1])
        B.append(RGB[2])

    cmap_dict = {}
    cmap_dict['red'] = [(x0[i], R[i], R[i]) for i in range(len(R))]
    cmap_dict['green'] = [(x0[i], G[i], G[i]) for i in range(len(G))]
    cmap_dict['blue'] = [(x0[i], B[i], B[i]) for i in range(len(B))]
    mymap = LinearSegmentedColormap('mymap', cmap_dict)
    return mymap


# --- polygon box
def draw_bbox(ax, bb):
    # boxstyle=square with pad=0, i.e. bbox itself.
    p_bbox = FancyBboxPatch((bb.xmin, bb.ymin),
                            abs(bb.width), abs(bb.height),
                            boxstyle="square,pad=0.",
                            linewidth=0,
                            edgecolor=None, facecolor="#DDDDDD", alpha=0.5,
                            zorder=1)  # zorder=1 as the background first
    ax.add_patch(p_bbox)


def set_xaxis_ticks(ax, xticks,direc='in'):
    ax.xaxis.set_major_locator(plt.MultipleLocator(xticks[0]))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(xticks[1]))
    ax.xaxis.set_tick_params(which='both', direction=direc)


def set_yaxis_ticks(ax, yticks,direc='in'):
    ax.yaxis.set_major_locator(plt.MultipleLocator(yticks[0]))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(yticks[1]))
    ax.yaxis.set_tick_params(which='both', direction=direc)
    ax.yaxis.set_ticks_position('both')


def freqcontour(ax, x, y, z, cmap='afmhot_r', vmin=-3, vmax=0, cb='on'):
    """
    plot log10-scale freq. contour  over a range of x, i.e. Ustar.

    inputs:
    ax - the figure axis
    x - array of parameter x (i.e. Ustar)
    y - array of parameter y (i.e. freq)
    z - array of parameter p (i.e. power)
    vmin - min of p range (-3 in default for log10 scale)
    vmax - max of p range (0 in default for log10 scale)
    cb - toggle 'on' or 'off' of the colorbar

    """
    N = 400
    for k in range(0,np.size(z)):
        if (z[k] == 0):
            z[k] = 1.0e-10
    xi = np.linspace(x.min(), x.max(), N)
    yi = np.linspace(y.min(), y.max(), N)
    zi = scipy.interpolate.griddata((x, y),
                                    z, (xi[None, :], yi[:, None]),
                                    method='linear')
    for k in range(N):
        zi[:,k] = (zi[:,k] / np.max(zi[:,k]))
    zi = np.log10(zi)
    ax.contourf(xi, yi, zi, 50, cmap=cmap,vmin=vmin, vmax=vmax)
#    ax.imshow(zi, vmin=vmin, vmax=vmax, cmap=cmap, origin="lower",
#              aspect="auto", extent=[x.min(), x.max(), y.min(),
#                                     y.max()])
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    if cb == 'on' or cb == 'On':
        axins = inset_axes(ax,
                           width="2.5%",  # width = 10% of parent_bbox width
                           height="100%",  # height : 50%
                           loc=6,
                           bbox_to_anchor=(1.02, 0., 1, 1),
                           bbox_transform=ax.transAxes,
                           borderpad=0)
        # get rid of unwanted colorbar range!
        cbar = matplotlib.colorbar.ColorbarBase(axins, cmap=cmap, norm=norm)
        axins.get_yaxis().set_tick_params(direction='in',
                                          which='both', width=0.5,
                                          labelleft='off')

    elif cb == 'off' or cb == 'Off':
        pass
    else:
        print('Toggle the colorbar \'on\' or \'off\'.')
        
def freqcontour1(ax, x, y, z, cmap='afmhot_r', vmin=-3, vmax=0, cb='on'):
#def freqcontour1(ax, x, y, z, cmap='gray_r', vmin=-3, vmax=0, cb='on'):
    """
    plot log10-scale freq. contour  over a range of x, i.e. Ustar.

    inputs:
    ax - the figure axis
    x - array of parameter x (i.e. Ustar)
    y - array of parameter y (i.e. freq)
    z - array of parameter p (i.e. power)
    vmin - min of p range (-3 in default for log10 scale)
    vmax - max of p range (0 in default for log10 scale)
    cb - toggle 'on' or 'off' of the colorbar

    """

    for k in range(0,np.size(z)):
        if (z[k] == 0):
            z[k] = 1.0e-4
    z = np.log10(z)
    z = np.reshape(z,(48,401)).T
    x = x[400:19248:401]
    y = y[0:401]
    print(np.shape(z),np.shape(x),np.shape(y))
    
#    ax.tricontourf(x, y, z, 20, cmap=cmap,vmin=vmin, vmax=vmax)
    ax.contourf(x, y, z, 50, cmap=cmap,vmin=vmin, vmax=vmax)
#    ax.imshow(z, vmin=vmin, vmax=vmax, cmap=cmap, origin="lower",
#              aspect="auto", extent=[x.min(), x.max(), y.min(),
#                                     y.max()])
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    if cb == 'on' or cb == 'On':
        axins = inset_axes(ax,
                           width="2.5%",  # width = 10% of parent_bbox width
                           height="100%",  # height : 50%
                           loc=6,
                           bbox_to_anchor=(1.02, 0., 1, 1),
                           bbox_transform=ax.transAxes,
                           borderpad=0)
        # get rid of unwanted colorbar range!
        cbar = matplotlib.colorbar.ColorbarBase(axins, cmap=cmap, norm=norm)
        axins.get_yaxis().set_tick_params(direction='in',
                                          which='both', width=0.5,
                                          labelleft='off')

    elif cb == 'off' or cb == 'Off':
        pass
    else:
        print('Toggle the colorbar \'on\' or \'off\'.')
        
def fft(signal,length,sampling_Freq):
    
    DFT = np.fft.rfft(signal,length)/length
    amp = 2*abs(DFT[0:int(length/2+1)]);
    
    freq = sampling_Freq/2*np.linspace(0,1,int(length/2+1));
    
    return amp,freq
    
def crossings_nonzero_pos2neg(data):
    pos = data > 0
    return (pos[:-1] & ~pos[1:]).nonzero()[0]
    
def crossings_nonzero_neg2pos(data):
    pos = data < 0
    return (pos[:-1] & ~pos[1:]).nonzero()[0]

def zero_cross(x, y, direction):
    if direction == 1:
        I = crossings_nonzero_pos2neg(y)
    else:
        I = crossings_nonzero_neg2pos(y)
    x1 = x[I]; y1 = abs(y[I]); x2 = x[I+1]; y2 = abs(y[I+1]);
    
    return (x1*y2+x2*y1)/(y1+y2)

def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)

def tecplot_reader1(file, nb_var):
    """Tecplot reader."""
    pos1 = 0
    nc = 0
    with open(file, 'r') as a:
        for idx, line in enumerate(a.readlines(),-3):
            if idx < 0:
                if line[1:5] == "ZONE":
                    N = int(re.findall(r',N=.*,', line)[0].strip(',N='))
                    arrays = np.zeros((N,nb_var), dtype=float)
            else:
                if (pos1 == N):
                    nc = nc+1
                    pos1 = 0
                    continue
                if (nc == nb_var):
                    break
                temp = np.fromstring(line, sep=' ')
                pos2 = len(temp) + pos1
                arrays[pos1:pos2,nc] = temp
                pos1 = pos2
    return arrays

# The following rotate function performs a rotation of the point point by the
# angle angle (counterclockwise, in radians) around origin, in the Cartesian 
# plane, with the usual axis conventions: x increasing from left to right, y 
# increasing vertically upwards. All points are represented as length-2 tuples
# of the form (x_coord, y_coord)
def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


if __name__=="__main__":
    from matplotlib.pyplot import plot, scatter, show
    series = [0,0,0,2,0,0,0,-2,0,0,0,2,0,0,0,-2,0]
    maxtab, mintab = peakdet(series,.3)
    plot(series)
    scatter(array(maxtab)[:,0], array(maxtab)[:,1], color='blue')
    scatter(array(mintab)[:,0], array(mintab)[:,1], color='red')
    show()
    
