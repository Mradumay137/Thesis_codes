#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 17:41:27 2023

@author: admin
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
from astropy.convolution import convolve, Box1DKernel
from scipy import integrate
from scipy.integrate import quad

for j in range(1,2):
    if j==1:
        s=str('S6.dat')
    if j==2:
        s=str('S6olddata1.dat')
    
    file=np.loadtxt(s,usecols=range(1,10),skiprows=2)
    
    col1=file[:,4]
    col2=file[:,5]
    col3=file[:,6]
    colav1=(col2+col3)/2
    
    col4=[]
    col5=[]
    
    xwoline=[]
    ywoline=[]
    
    def func2(x, a, b, c,d,e):
        return a*x**4+b*x**3+c*x**2+d*x+e
    
    xwline=[]
    ywline=[]
    for i in range(len(col1)):
        if col1[i]>1667 and col1[i]<1668:
            xwline.append(col1[i])
            ywline.append(colav1[i])
        
        if col1[i]>1667 and col1[i]<1668:
            if col1[i]>1667.38 and col1[i]<1667.44:
                continue
            if col1[i]>1667.73 and col1[i]<1667.80:
                continue
            ywoline.append(colav1[i])
            xwoline.append(col1[i])
    
    popt, pcov = curve_fit(func2, xwoline, ywoline)
    
    ywline=np.asarray(ywline)
    xwline=np.asarray(xwline)
    
    ysub=ywline-func2(xwline, *popt)
        
    
    for i in range(len(xwline)):
        if xwline[i]> 1667.2 and xwline[i]<1667.6:
            col4.append(xwline[i])
            col5.append(ysub[i])
        
    def func(x,a,b,c,d,e,f):
        return a*np.exp(-(x-b)**2/(2*c**2))+d*np.exp(-(x-e)**2/(2*f**2))
    velocity=[]
    c=3*10**5
    for i in range(0,len(col4)):
        v=((1667.3590-col4[i])*c/col4[i])
        velocity.append(v)
    
    col5smooth=convolve(col5,Box1DKernel(7))
    velocity=np.asarray(velocity)
    
    popt, pcov = curve_fit(func, velocity,col5smooth,p0=(0.04,-7,0.1,0.05,-7,0.2),method='lm')
    
    a=popt[0]
    b=popt[1]
    c=popt[2]
    d=popt[3]
    e=popt[4]
    f=popt[5]
    
    NOH=quad(func,velocity[-1],velocity[0],args=(a,b,c,d,e,f))
    print(NOH)
    col6=col5smooth[np.where((velocity>-20)&(velocity<5))]
    velocitychop=velocity[np.where((velocity>-20)&(velocity<5))]
    
    NOH2=integrate.simpson(col6,dx=velocitychop[0]-velocitychop[1])
    print(NOH2)
    
    plt.plot(velocity, func(velocity,*popt),label=s)
    
    plt.plot(velocity,col5smooth,label=s)
    
    plt.legend()
    
