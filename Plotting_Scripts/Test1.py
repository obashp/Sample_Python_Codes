import os
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import Atul_myplt as myplt
import Utils as ut

folder_path = 'E:\BACKUP_NIHAR_1_OCT_2020\\Nihar\\VIV_tests\\D_O\\'
Ur = 9.5

x1, y1 = ut.phasedata('2y', folder_path, 'D_O_', Ur, 0, 5)
plt.plot(x1[0:1500],y1[0:1500],'k',markersize=1.0)
y2 = y1
x1, y1 = ut.phasedata('2y', folder_path, 'D_O_', Ur, 0, 4)
# plt.plot(x1[0:1500],y1[0:1500],'b',markersize=1.0)
x1, y1 = ut.phasedata('2y', folder_path, 'D_O_', Ur, 0, 3)
plt.plot(x1[0:1500],y2[0:1500]-0.25*3.14159*y1[0:1500],'c',markersize=1.0)
x1, y1 = ut.phasedata('2y', folder_path, 'D_O_', Ur, 0, 1)
plt.plot(x1[0:1500],y1[0:1500],'g',markersize=1.0)