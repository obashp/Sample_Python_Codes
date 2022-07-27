# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 11:52:47 2020

@author: Laxman
"""


import os
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import gridspec
import matplotlib.lines as mlines


folder_path0 = 'F:\\Nihar\\VIV_tests\\D_O'
Ampfilelist = os.path.join(folder_path0,'Afilelist.dat')
with open(Ampfilelist,'r') as f:
    filelist = np.loadtxt(f,dtype=str)

nfiles = len(filelist) - 1
nUr = float(filelist[nfiles])    
print(nfiles,nUr)
