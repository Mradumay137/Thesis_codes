#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:07:45 2023

@author: admin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

file=np.loadtxt("0056-0717_1446.txt",usecols=range(1,3))
file1=np.loadtxt("H1_S6.dat",usecols=range(1,10))

col1=file1[:,4]
col2=file1[:,5]
col3=file1[:,6]
colav1=(col2+col3)/2

col4=file[:,0]
col5=file[:,1]

reducedx=[]
reducedy=[]

for i in range(len(col1)):
    if col1[i]>1415 and col1[i]<1425:
        if col1[i]>1420 and col1[i]<1421:
            continue
        reducedx.append(col1[i])
        reducedy.append(colav1[i])

reducedx1=[]        
reducedy1=[]
for i in range(len(col1)):
    if col1[i]>1415 and col1[i]<1425:
        reducedx1.append(col1[i])
        reducedy1.append(colav1[i])

reducedx1=np.asarray(reducedx1)
def func(x, a, b, c,d,e,f):
    return a*(np.sin(x))**5+b*(np.sin(x))**4+c*(np.sin(x))**3+d*(np.sin(x))**2+e*np.sin(x)+f

popt, pcov = curve_fit(func, reducedx, reducedy)
subtracted=reducedy1-func(reducedx1, *popt)

velocity=[]
c=3*10**5
for i in range(0,len(reducedx1)):
    v=((1420.406-reducedx1[i])*c/reducedx1[i])
    velocity.append(v)      
    
plt.plot(col4,col5,label='GASS data')
plt.plot(velocity, subtracted,label='S6 data')
plt.legend()
plt.title('H1_comparison')