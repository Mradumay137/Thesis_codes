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
from scipy.integrate import quad
from scipy import integrate

plt.figure(figsize=(15,15))
for n in range(3,5):
    
    if n==1:
        s=str('CHmet1.dat')
    if n==2:
        s=str('CHmethod2.dat')  
    if n==3:
        
        s=str('S6CH.dat')
    if n==4:
        s=str('S6.dat')   
        
    if n==1:
        p=str('ON-OFF/OFF (Method 1)')
   
    if n==2:
       p=str('ON-OFF with noise cal (Method 2)')  
    
    if n==3:
       p=str('ON with 1934')
    
    if n==4:
        p=str('ON with 1934 (Method 4)')    
        
    
    file1=np.loadtxt(s,usecols=range(1,10))

    col1=file1[:,4]
    col2=file1[:,5]
    col3=file1[:,6]
    colav1=(col2+col3)/2
    
    if n==3:
        colav1=colav1*0.26
    
    yvalue1=[]
    xvalue1=[]
    yvalue4=colav1[np.where((col1>3335)&(col1<3336))]
    xvalue4=col1[np.where((col1>3335)&(col1<3336))]
    for i in range(len(col1)):
        if col1[i]>3335 and col1[i]<3336:
            if col1[i]>3335.5 and col1[i]<3335.62:
                continue
            yvalue1.append(colav1[i])
            xvalue1.append(col1[i])
    
    
    def func(x, a, b, c,d,e,f,g,h,i,j,k,l):
        return a*x**11+b*x**10+c*x**9+d*x**8+e*x**7+f*x**6+g*x**5+h*x**4+i*x**3+j*x**2 + k*x + l
    
    
    def func1(x, a, b,c,d,e):
        return a*x**3+b*x**2+d*x+e
    def func4(x, a, b,c,d,e,f):
        return a*x**4+b*x**3+d*x**2+e*x+f
   
    def func2(x, a, b, c,d,e,f,g):
        return a*x**6+b*x**5+c*x**4+d*x**3+e*x**2+f*x+g
    
    def polyN(x, p):
    	y = 0
    	i = len(p) - 1
    	p = np.flip(p)
    	while (i >= 0):
    		#print(i)
    		y += p[i]*x**float(i)
    		i -= 1
    	return y
    xvalue1=np.asarray(xvalue1)
    yvalue1=np.asarray(yvalue1) 
    p = np.polyfit(xvalue1, yvalue1, 6)

    y_model = polyN(xvalue4, p)
    
    
    
    """
    popt, pcov = curve_fit(polyN, xvalue1, yvalue1,maxfev=1000)
    """
    xvalue3=col1[np.where((col1>3335.51)&(col1<3335.61))]
    yvalue3=colav1[np.where((col1>3335.51)&(col1<3335.61))]
    
    xvalue2=col1[np.where((col1>3335)&(col1<3336))]
    yvalue2=colav1[np.where((col1>3335)&(col1<3336))]
    
    
    convol1=convolve(yvalue2,Box1DKernel(1))
    if n==3:
        convol1=convolve(yvalue2,Box1DKernel(2))
    
    convol1chopped=convol1[10:len(convol1)-10]
    xvalue2chop=xvalue2[10:len(xvalue2)-10]
    
    
   
    """
    yvaluesub=polyN(xvalue1, 4) 
    """
    #plt.plot(xvalue4,yvalue4,label='with')
    #plt.plot(xvalue1,yvalue1,label='wl') 
    #plt.plot(xvalue4,y_model,label='model')
    
    
    velocity=[]
    c=3*10**5
    for i in range(0,len(xvalue2)):
        v=((3335.4810-xvalue2[i])*c/xvalue2[i])
        velocity.append(v)
           
    velocity=np.asarray(velocity) 
    
    NCH2=integrate.simpson(convol1-y_model,dx=velocity[0]-velocity[1])
    
    print(NCH2)
    plt.title("3335 MHz spectral line",fontsize=20)
    
    
    plt.plot(velocity,yvalue2-y_model,label='ON with 1934 (Method 4)')
    plt.axhline(y=0,color='k',label='Zero Flux line')
    plt.xlabel("Frequency (MHz)",fontsize=20)
    plt.ylabel("Flux (Jansky)",fontsize=20)
    plt.legend(fontsize=20)
    plt.xlim(-40,30)
    
    plt.savefig('CHmet4',dpi=300)
    
    print(np.std(convol1-y_model))
    
    """
    mask = ((xvalue2chop > 3335.4) & (xvalue2chop < 3335.51)) | ((xvalue2chop > 3335.61) & (xvalue2chop < 3335.7))

    yvaluesubsampl = yvaluesub[mask]
    
    
    popt4, pcov4 = curve_fit(func1, xvalue2chop,convol1chopped,maxfev=5000)
    
    yvaluesubwo=yvaluesub-func1(xvalue2chop,*popt4)
    
    yvalue4=colav1[np.where((col1>3335)&(col1<3335.3))]
    xvalue4=col1[np.where((col1>3335)&(col1<3335.3))]
    
    yvaluesub1=yvalue4-func(xvalue4, *popt)
    
    
    def func4(x,a,b,c,d,e,f):
        return a*np.exp(-(x-b)**2/(2*c**2))+d*np.exp(-(x-e)**2/(2*f**2))
        
   
    
      
    
    
    popt2, pcov2 = curve_fit(func4, velocity,yvaluesub,p0=(0.04,-7,0.2,0.05,-8,0.2),method='lm')
     
    a=popt2[0]
    b=popt2[1]
    c=popt2[2]
    d=popt2[3]
    e=popt2[4]
    f=popt2[5]
     
    NCH=quad(func4,velocity[-1],velocity[0],args=(a,b,c,d,e,f))
    
        
    Nextra=integrate.simpson(yvaluesubsampl,dx=velocity[1]-velocity[0])
    
    
    print(NCH2)
    
    


    plt.figure(figsize=(7,15))
    
    plt.plot(xvalue2,yvalue2,label="Method {0}".format(n))
    

    plt.subplot(3,1,1)
    plt.plot(xvalue2,yvalue2,'k',label='On-source spectrum before baseline correction')
    plt.subplot(3,1,2)
    
    
    
    if n==4 or n==3:
        plt.axhline(y=0,color='k',label='Zero Flux line')
    
    plt.xlim(-40,30)
    
    
    yvaluesubhist=yvalue3-func2(xvalue3, *popt)
    
    suma=0
    
    for beta in range(0,len(yvaluesubhist)):
        
        suma=suma+(yvaluesubhist[beta])**2
        
    
    alpha=integrate.simpson(yvaluesubhist,dx=velocity[0]-velocity[1])
    print(alpha)
    mask2 = ((xvalue2chop > 3335.44) & (xvalue2chop < 3335.66)) 

    yvalueplot = yvaluesub[mask2]
    xvalueplot = velocity[mask2]
    


plt.xlabel("Frequency (MHz)",fontsize=20)
plt.ylabel("Flux (Jansky)",fontsize=20)
plt.legend(fontsize=20)

hist,binedge=np.histogram(yvaluesubhist,bins=10,density=False)
weights=np.ones_like(yvaluesubhist)/len(yvaluesubhist)

plt.hist(yvaluesubhist, weights=weights)

  
plt.title("3335 MHz spectral line",fontsize=20)
plt.savefig('CHmet4',dpi=300)
"""