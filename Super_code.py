#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 13:55:08 2023

@author: admin
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from scipy.optimize import curve_fit
n=1
plt.figure(figsize=(30,20))
for m in range(1,7):
    if n==1:
        s=str('2a_scal.dat')
    if n==2:
        s=str('2b_tcal.dat')  
    if n==3:
        s=str('2c_scal.dat')
    if n==4:
        s=str('2c_tcal.dat')
    if n==5:
        s=str('2d_scal.dat')
    if n==6: 
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
    yvalueoriginal2=[]
    xvalueoriginal2=[]
    yvalueoriginal3=[]
    xvalueoriginal3=[]
    
    for i in range(len(col1)):
        
        if col1[i]>1665.2 and col1[i]<1665.7:
            yvalueoriginal2.append(colav1[i]) 
            xvalueoriginal2.append(col1[i])
        if col1[i]>1667 and col1[i]<1668:
            yvalueoriginal3.append(colav1[i]) 
            xvalueoriginal3.append(col1[i])
           
        
        if col1[i]>1612.1 and col1[i]<1612.5:
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
        
        if col1[i]>1720.1 and col1[i]<1721:
            yvalue4.append(colav1[i]) 
            xvalue4.append(col1[i])
     
    velocity1=[]
    velocity2=[]
    velocity3=[]
    velocity4=[]
    
    c=3*10**5
    
    for i in range(0,len(xvalue1)):
        v1=((1612.2310-xvalue1[i])*c/xvalue1[i])
        velocity1.append(v1)
    for i in range(0,len(xvalueoriginal2)):     
        v2=((1665.4019-xvalueoriginal2[i])*c/xvalueoriginal2[i])
        velocity2.append(v2)
    for i in range(0,len(xvalueoriginal3)):
        v3=((1667.3590-xvalueoriginal3[i])*c/xvalueoriginal3[i])
        velocity3.append(v3)
    for i in range(0,len(xvalue4)):
        v4=((1720.5300-xvalue4[i])*c/xvalue4[i])
        velocity4.append(v4)
    
    convol1=convolve(yvalue1,Box1DKernel(9))
    convol1chopped=convol1[10:len(convol1)-10]
    
    convol2=convolve(yvalueoriginal2,Box1DKernel(9))
    convol2chopped=convol2[10:len(convol2)-10]
    
    convol3=convolve(yvalueoriginal3,Box1DKernel(9))
    convol3chopped=convol3[10:len(convol3)-10]
    
    convol4=convolve(yvalue4,Box1DKernel(9))
    convol4chopped=convol4[10:len(convol4)-10]
    
    def func2(x, a, b, c,d):
        return a*x**3+b*x**2+c*x+d
    
    def func(x, a, b, c,d,e):
        return a*x**4+b*x**3+c*x**2+d*x**1+e
    
    xvalue1=np.asarray(xvalue1)
    xvalue2=np.asarray(xvalue2)
    xvalue3=np.asarray(xvalue3)
    xvalue4=np.asarray(xvalue4)
    xvalueoriginal2=np.asarray(xvalueoriginal2)
    xvalueoriginal3=np.asarray(xvalueoriginal3)
   
    
   
    xvalueoriginalchop2=xvalueoriginal2[10:len(xvalueoriginal2)-10] 
   
    if n==5 or n==6:
        popt2, pcov2 = curve_fit(func2, xvalue2, yvalue2)
        
        yvaluesub2=convol2chopped-func2(xvalueoriginalchop2, *popt2)
    else:
        popt2, pcov2 = curve_fit(func, xvalue2, yvalue2)
        yvaluesub2=convol2chopped-func(xvalueoriginalchop2, *popt2)
    
    if n==5 or n==6:
        popt1, pcov1 = curve_fit(func2, xvalue1, yvalue1)
        xvalue1=xvalue1[10:len(xvalue1)-10]
        yvaluesub1=convol1chopped-func2(xvalue1, *popt1)
    else:
        popt1, pcov1 = curve_fit(func, xvalue1, yvalue1)
        xvalue1=xvalue1[10:len(xvalue1)-10]
        yvaluesub1=convol1chopped-func(xvalue1, *popt1)
    
    
    popt3, pcov3 = curve_fit(func, xvalue3, yvalue3)
    xvalueoriginalchop3=xvalueoriginal3[10:len(xvalueoriginal3)-10]
    yvaluesub3=convol3chopped-func(xvalueoriginalchop3, *popt3)
    
    popt4, pcov4 = curve_fit(func, xvalue4, yvalue4)
    xvalue4=xvalue4[10:len(xvalue4)-10]
    yvaluesub4=convol4chopped-func(xvalue4, *popt4)
    
    velocitychop1=velocity1[10:len(velocity1)-10]
    velocitychop1=np.asarray(velocitychop1)
    velocitychop2=velocity2[10:len(velocity2)-10]
    velocitychop2=np.asarray(velocitychop2)
    velocitychop3=velocity3[10:len(velocity3)-10]
    velocitychop3=np.asarray(velocitychop3)
    velocitychop4=velocity4[10:len(velocity4)-10]
    velocitychop4=np.asarray(velocitychop4)
    
        
    plt.subplot(6, 1, n)    
    plt.plot(velocitychop1,yvaluesub1,label='1612')
    plt.plot(velocitychop2,yvaluesub2,label='1665')
    plt.plot(velocitychop3,yvaluesub3,label='1667')
    plt.plot(velocitychop4,yvaluesub4,label='1720')
    plt.axhline(y=0,color='black',linestyle='dashdot')
    plt.xlabel("Velocity(km/s)")
    plt.ylabel("Flux")
    plt.xlim(-20,5)
    plt.ylim(-0.05,0.1)
    
    """
    plt.subplot(4, 1, 2)
    plt.plot(velocitychop2,yvaluesub2,label=s)
    plt.axhline(y=0,color='black',linestyle='dashdot')
    plt.xlabel("Velocity(km/s)")
    plt.ylabel("Flux")
    plt.xlim(-20,5)
    plt.ylim(-0.1,0.1)
    
    plt.subplot(4, 1, 3)
    plt.plot(velocitychop3,yvaluesub3,label=s)
    plt.axhline(y=0,color='black',linestyle='dashdot')
    plt.xlabel("Velocity(km/s)")
    plt.ylabel("Flux")
    plt.xlim(-20,5)
    plt.ylim(-0.1,0.1)
    
    plt.subplot(4, 1, 4)
    plt.plot(velocitychop4,yvaluesub4,label=s)
    plt.axhline(y=0,color='black',linestyle='dashdot')
    plt.xlabel("Velocity(km/s)")
    plt.ylabel("Flux")
    plt.xlim(-20,5)
    plt.ylim(-0.1,0.1)
    """
    n=n+1

plt.legend()
    