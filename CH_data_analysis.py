#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 17:51:32 2023

@author: admin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.convolution import convolve, Box1DKernel

for n in range(1,4):
    
    if n==1:
        s=str('2a_scal_CH_S6.dat')
    if n==2:
        s=str('2c_scal_CH_S6.dat')  
    if n==3:
        s=str('2d_scal_CH_S6.dat')
    
    file1=np.loadtxt(s,usecols=range(1,10))

    col1=file1[:,4]
    col2=file1[:,5]
    col3=file1[:,6]
    colav1=(col2+col3)/2
    
    yvalue1=[]
    xvalue1=[]
    
    for i in range(len(col1)):
        if col1[i]>3335 and col1[i]<3336:
            if col1[i]>3335.4 and col1[i]<3335.7:
                continue
            yvalue1.append(colav1[i])
            xvalue1.append(col1[i])
    
    def func(x, a, b, c,d,e):
        return a*x**4+b*x**3+c*x**2+d*x**1+e
   
    def func2(x, a, b, c,d,e,f,g):
        return a*x**6+b*x**5+c*x**4+d*x**3+e*x**2+f*x+g
        
    popt, pcov = curve_fit(func, xvalue1, yvalue1)
    
    xvalue2=col1[np.where((col1>3335)&(col1<3336))]
    yvalue2=colav1[np.where((col1>3335)&(col1<3336))]
    
    convol1=convolve(yvalue2,Box1DKernel(5))
    
    convol1chopped=convol1[10:len(convol1)-10]
    xvalue2chop=xvalue2[10:len(xvalue2)-10]
    
    yvaluesub=convol1chopped-func(xvalue2chop, *popt)
        
    plt.plot(xvalue2chop,yvaluesub,label="Method {0}".format(n))
    plt.legend()