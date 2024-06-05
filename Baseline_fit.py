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

plt.figure(figsize=(10,10))
for m in range(0,6):
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
    
    file1=np.loadtxt(s,usecols=range(1,10))
    
    col1=file1[:,4]
    col2=file1[:,5]
    col3=file1[:,6]
    colav1=(col2+col3)/2
    yvalue1=[]
    yvalue2=[]
    yvalue3=[]
    yvalue4=[]
    xvalue1=[]
    xvalue2=[]
    xvalue3=[]
    xvalue4=[]
    yvalueoriginal=[]
    xvalueoriginal=[]
    
    """
    xvalue1= col1[np.where((col1>1610) & (col1<1614))]
    yvalue1=colav1[np.where((col1>1610) & (col1<1614))]
    xvalue2= col1[np.where((col1>1665) & (col1<1666))]
    yvalue2=colav1[np.where((col1>1665) & (col1<1666))]
    xvalue3= col1[np.where((col1>1667) & (col1<1668))]
    yvalue3=colav1[np.where((col1>1667) & (col1<1668))]
    """
    
    for i in range(len(col1)):
        if col1[i]>1665.2 and col1[i]<1665.7:
            yvalueoriginal.append(colav1[i]) 
            xvalueoriginal.append(col1[i])
        if col1[i]>1610 and col1[i]<1614:
            yvalue1.append(colav1[i])
            xvalue1.append(col1[i])
        if col1[i]>1665.2 and col1[i]<1665.7:
            if col1[i]>1665.42 and col1[i]<1665.46:
                continue
            yvalue2.append(colav1[i])
            xvalue2.append(col1[i])
            
        
           
        if col1[i]>1667 and col1[i]<1668:
            if col1[i]>1667.38 and col1[i]<1667.44:
                continue
            if col1[i]>1667.73 and col1[i]<1667.80:
                continue
            yvalue3.append(colav1[i])
            xvalue3.append(col1[i])
        
        
        if col1[i]>1718 and col1[i]<1721:
            yvalue4.append(colav1[i]) 
            xvalue4.append(col1[i])
        
    def func(x, a, b, c,d,e,f):
        return a*x**5+b*x**4+c*x**3+d*x**2+e*x**1+f
    
    def func2(x, a, b, c,d,e):
        return a*x**4+b*x**3+c*x**2+d*x**1+e
    
    xvalue2=np.asarray(xvalue2)
    xvalueoriginal=np.asarray(xvalueoriginal)
    xvalue2chop=xvalue2[10:len(xvalue3)-10]
    yvalue3chop=yvalue3[10:len(xvalue3)-10]
    
    convol1=convolve(yvalue3,Box1DKernel(7))
    
    if n==4 or n==5:
        popt, pcov = curve_fit(func, xvalue2, yvalue2)
        
    else:
        popt, pcov = curve_fit(func2, xvalue2, yvalue2)
    convol1=convolve(yvalueoriginal,Box1DKernel(7))
    convol1chopped=convol1[10:len(convol1)-10]
    xvalueoriginalchop=xvalueoriginal[10:len(xvalueoriginal)-10]
    if n==4 or n==5:
        yvaluesub=convol1chopped-func(xvalueoriginalchop, *popt)
    else:
        yvaluesub=convol1chopped-func2(xvalueoriginalchop, *popt)
    """
    plt.plot(xvalue2, yvalue2)
    plt.plot(xvalue2, func(xvalue2, *popt))
    """
    
    """
    plt.plot(xvalue2, func(xvalue2, *popt))
    plt.plot(xvalue3chop,convol1chopped)
    """
    
    velocity=[]
    c=3*10**5
    for i in range(0,len(xvalueoriginal)):
        v=((1665.4019-xvalueoriginal[i])*c/xvalueoriginal[i])
        velocity.append(v)
    velocitychop=velocity[10:len(velocity)-10]
    velocitychop=np.asarray(velocitychop)
    if n==0 or n==2 or n==4:
        
        plt.subplot(2, 1, 1)    
        plt.plot(velocitychop,yvaluesub,label=s)
        plt.axhline(y=0,color='black',linestyle='dashdot')
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Flux")
        plt.xlim(-20,5)
        plt.title("Calibration methods")
        plt.legend()
    
    if n==1 or n==3 or n==5:
        
        plt.subplot(2, 1, 2)
        plt.plot(velocitychop,yvaluesub,label=s)
        plt.axhline(y=0,color='black',linestyle='dashdot')
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Flux")
        plt.xlim(-20,5)
        
        plt.legend()
    n=n+1
"""
print(convol1[0:10])

print(convol1Chopped[0:10])

col4=file2[:,4]
col5=file2[:,5]
col6=file2[:,6]
colav2=(col5+col6)/2

col7=file3[:,4]
col8=file3[:,5]
col9=file3[:,6]
colav3=(col8+col9)/2

col10=file4[:,4]
col11=file4[:,5]
col12=file4[:,6]
colav4=(col11+col12)/2

col13=file5[:,4]
col14=file5[:,5]
col15=file5[:,6]
colav5=(col14+col15)/2

plt.plot(col1,(col2+col3)/2,'y',label='2a_scal')
plt.plot(col4,(col5+col6)/2,'r',label='2c_scal')
plt.plot(col7,(col8+col9)/2,'b',label='2c_tcal')
plt.plot(col10,(col11+col12)/2,'g',label='2d_scal')
plt.plot(col13,(col14+col15)/2,'k',label='2d_tcal')
"""