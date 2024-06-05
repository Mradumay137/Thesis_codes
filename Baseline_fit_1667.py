#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 14:37:38 2023

@author: admin
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 16:09:35 2023

@author: admin
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.convolution import convolve, Box1DKernel
from scipy.optimize import curve_fit
import pandas as pd

n=0

file1=np.loadtxt("2a_scal.dat",usecols=range(1,10))
file2=np.loadtxt("2c_scal.dat",usecols=range(1,10))
file3=np.loadtxt("2c_tcal.dat",usecols=range(1,10))
file4=np.loadtxt("2d_scal.dat",usecols=range(1,10))
file5=np.loadtxt("2d_tcal.dat",usecols=range(1,10))
file6=np.loadtxt("2b_tcal.dat",usecols=range(1,10))
file7=np.loadtxt("2e_scal.dat",usecols=range(1,10))

plt.figure(figsize=(8,5))
for m in range(0,7):
    if n==0:
        s=str('2a_scal.dat')
    if n==1:
        s=str('2b_tcal.dat')  
    if n==2:
        s=str('2c_scal.dat')
    if n==3:
        s=str('2c_tcal.dat')
    if n==4:
        s=str('2d_scal.dat')
    if n==5: 
        s=str('2d_tcal.dat')
    if n==6:
        s=str('2e_scal.dat')
    
    file1=np.loadtxt(s,usecols=range(1,10))
    
    col1=file1[:,4]
    col2=file1[:,5]
    col3=file1[:,6]
    colav1=(col2+col3)/2
    yvalue3=[]
    xvalue3=[]
    xvalueoriginal=col1[np.where((col1>1667)&(col1<1668))]
    yvalueoriginal=colav1[np.where((col1>1667)&(col1<1668))]
    
    for i in range(len(col1)):
       if col1[i]>1667 and col1[i]<1668:
            if col1[i]>1667.38 and col1[i]<1667.44:
                continue
            if col1[i]>1667.73 and col1[i]<1667.80:
                continue
            yvalue3.append(colav1[i])
            xvalue3.append(col1[i])
        
    def func(x, a, b, c,d,e):
        return a*x**4+b*x**3+c*x**2+d*x**1+e
    
    xvalue3=np.asarray(xvalue3)
    xvalueoriginal=np.asarray(xvalueoriginal)
    
    popt, pcov = curve_fit(func, xvalue3, yvalue3)
        
    convol1=convolve(yvalueoriginal,Box1DKernel(7))
    convol1chopped=convol1[10:len(convol1)-10]
    xvalueoriginalchop=xvalueoriginal[10:len(xvalueoriginal)-10]
    
    yvaluesub=convol1chopped-func(xvalueoriginalchop, *popt)
    
    velocity=[]
    c=3*10**5
    for i in range(0,len(xvalueoriginal)):
        v=((1667.3590-xvalueoriginal[i])*c/xvalueoriginal[i])
        velocity.append(v)
    velocitychop=velocity[10:len(velocity)-10]
    velocitychop=np.asarray(velocitychop)
    if n==0 or n==2 or n==4 or n==6:
        """
        plt.subplot(2, 1, 1)  
        """
        plt.plot(velocitychop,yvaluesub*0.7,label=s)
        plt.axhline(y=0,color='black',linestyle='dashdot')
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Flux")
        plt.xlim(-20,5)
        plt.title("Calibration methods")
        plt.legend()
    
    if n==1 or n==3 or n==5:
        """
        plt.subplot(2, 1, 2)
        """
        plt.plot(velocitychop,yvaluesub,label=s)
        plt.axhline(y=0,color='black',linestyle='dashdot')
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Flux")
        plt.xlim(-20,5)
        
        plt.legend()
    n=n+1
    yvaluenew=yvaluesub[(np.where(velocitychop>0) and (velocitychop<5))]
    std=np.std(yvaluenew)
    print(std)

