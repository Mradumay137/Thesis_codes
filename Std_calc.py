#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 17:40:33 2024

@author: admin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

inttime=np.array([7461.7,7460.7,7251.8,10836.8,7254.8,7301.8,7463.7,7252.8,7250.8,7253.8,7249.8,7253.8,8092.3,8085.3,8079.3,8090.3,7853.5,7969.4,7872.5,7675.6,7766.5,7771.5,7671.6,7667.6,6508.2,6511.2,6510.2,6507.2,6505.3,6297.4,6302.4,6422.3,6627.2,8031.4,7070.9,6841.1,5980.8,7175.2,7175.2,7176.2])
std=[]

for i in range(1,2):
    """
    s=f'S{i}.dat'
    """
    s='stat4.dat'
    file1=np.loadtxt(s,usecols=range(1,11),skiprows=2)
    file2=np.loadtxt('S80H.dat',usecols=range(1,11),skiprows=2)
    
    col1=file1[:,4]
    col2=file1[:,5]
    col3=file1[:,6]
    colav1=(col2+col3)/2
    
    colav1red=colav1[np.where((col1>1667.5)& (col1<1668))]
    
    col4=file2[:,4]
    col5=file2[:,5]
    col6=file2[:,6]
    colav2=(col5+col6)/2
    
    colav2red=colav2[np.where((col4>1668)& (col4<1669))]
    
    m=stats.mode(colav1)
    
    sd1=np.std(colav1red)
    std.append(sd1)
    sd2=np.std(colav2red)

    print(sd1)
"""
plt.scatter(inttime,std )
plt.ylim(0,0.02)
plt.ylabel('rms noise')
plt.xlabel('Integration time(s)')
plt.savefig('RMs_fit.png', dpi=300)
"""