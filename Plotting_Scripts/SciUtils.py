import os,glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

def butter_lowpass(cutoff, fs, order=4):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def read_file(filename,idxlist):
    with open(filename,'r') as f:
        data = np.loadtxt(filename,skiprows=2)
        
    
    datalen = len(data[:,0])
    varis = np.zeros([datalen,len(idxlist)])
    print(datalen)
    
    i = 0
    while i < len(idxlist):
        varis[:,i] = data[:,idxlist[i]]
        i = i+1
    
    return varis

def overlap_idx(t1,t2,dt):
    stidx = 0
    i = 0
    while i < len(t2):
        if abs(t2[i] - t1[0]) < dt:
            stidx = i
            break
        i = i+1
    
    return stidx

def plot_ty(locy1y2,locy1cy1,locy2cy2,folder_path):
    viv_fname = 'viv_body2dir2_probe.dat'
    lift_fname = 'drag_lift_body_002.dat'
    
    vivfile = os.path.join(folder_path,viv_fname)
    liftfile = os.path.join(folder_path,lift_fname)
    
    vivvaris = read_file(vivfile,[1,2])
    liftvaris = read_file(liftfile,[0,6])
    
    edidx = len(vivvaris[:,0])-100
    stidx = int(0.5*edidx)
    span = edidx - stidx
    
    time = vivvaris[stidx:edidx,0]
    yO = vivvaris[stidx:edidx,1]
    
    dt = time[1]-time[0]
    fs = 1.0/dt
    
    lstidx = overlap_idx(time,liftvaris[:,0],dt)
    cyO = liftvaris[lstidx:lstidx+span,1]
    cyO = butter_lowpass_filtfilt(cyO,8.0,fs)
    
    viv_fname = 'viv_body1dir2_probe.dat'
    lift_fname = 'drag_lift_body_001.dat'
    
    vivfile = os.path.join(folder_path,viv_fname)
    liftfile = os.path.join(folder_path,lift_fname)
    
    vivvaris = read_file(vivfile,[1,2])
    tstidx = overlap_idx(time,vivvaris[:,0],dt)
    yD = vivvaris[tstidx:tstidx+span,1]
    
    liftvaris = read_file(liftfile,[0,6])
    lstidx = overlap_idx(time,liftvaris[:,0],dt)
    cyD = liftvaris[lstidx:lstidx+span,1]
    cyD = butter_lowpass_filtfilt(cyD,8.0,fs)
        
    ast = -7000
    aed = -3000
    
    t = time[ast:aed]
    yO = yO[ast:aed]
    cyO = cyO[ast:aed]
    yD = yD[ast:aed]
    cyD = cyD[ast:aed]
    
    locy1y2.plot(t,yO,c='xkcd:black',ls='solid',lw = 0.7)
    locy1y2.plot(t,yD,c='xkcd:blue',ls='solid',lw = 0.7)
    maxy = max(max(yO),max(yD))
    miny = min(min(yO),min(yD))
    y1y2lim = [miny-0.1,maxy+0.1]
    
    locy1cy1.plot(t,cyD,c='xkcd:black',ls='solid',lw=0.8)
    locy1cy1.plot(t,yD,c='xkcd:blue',ls='dashed',lw=0.8)
    maxyD = max(max(cyD),max(yD))
    minyD = min(min(cyD),min(yD))
    y1cy1lim = [minyD-0.15,maxyD+0.15]
    
    
    locy2cy2.plot(t,cyO,c='xkcd:black',ls='solid',lw=0.8)
    locy2cy2.plot(t,yO,c='xkcd:blue',ls='dashed',lw=0.8)
    maxyO = max(max(cyO),max(yO))
    minyO = min(min(cyO),min(yO))
    y2cy2lim = [minyO-0.15,maxyO+0.15]
    
    AyD = [min(yD),max(yD)]
    AyO = [min(yO),max(yO)]
    
    tlim = [t[0],t[-1]]
    return (tlim,y1y2lim,y1cy1lim,y2cy2lim,AyD,AyO)





        