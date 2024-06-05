#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 14:59:45 2023

@author: admin
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.convolution import convolve, Box1DKernel

file=np.loadtxt("1a.dat",usecols=range(1,10))
file1=np.loadtxt("1b.dat",usecols=range(1,10))

col1=file[:,4]
col2=file[:,5]
col3=file[:,6]

col4=file1[:,4]
col5=file1[:,5]
col6=file1[:,6]

col7=[]
col8=[]
 
for u in range(len(col1)):
    if col1[u]>1664.2 and col1[u]<1667:
        col7.append(col3[u])
        
        
for u in range(len(col4)):
    if col4[u]>1664.2 and col4[u]<1667:
        col8.append(col6[u])


std1a=np.std(col7)
std1b=np.std(col8)

print(std1a,std1b)

convol1=convolve(col3,Box1DKernel(13))
convol2=convolve(col6,Box1DKernel(13))

plt.figure(figsize=(10,6))

plt.plot(col1,convol1,'k',label='Approach 1')
plt.plot(col4,convol2,label='Approach 2')
plt.legend()
plt.xlabel("Frequency (MHz)",fontsize=15)
plt.ylabel("Temperature (K)",fontsize=15)

plt.xlim(1666,1668)
plt.ylim(-0.01,0)

plt.savefig('Order',dpi=300)
plt.show()
