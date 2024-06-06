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
import csv
I=[]
errorval=[]
for omicron in range(6,11):
    plt.figure(figsize=(15,15))
    for n in range(1,5):
        if omicron==9:
            continue
        
        if n==1:
            s=str(f'S{omicron}met1CH.dat')
        if n==2:
            s=str(f'S{omicron}met2CH.dat')  
        if n==3:
            
            s=str(f'S{omicron}met3CH.dat')
        if n==4:
            s=str(f'S{omicron}CH.dat')   
            
        if n==1:
            po=str('ON-OFF/OFF (Method 1)')
       
        if n==2:
           po=str('ON-OFF with noise cal (Method 2)')  
        
        if n==3:
           po=str('ON with noisecal (Method 3)')
        
        if n==4:
            po=str('ON with 1934 (Method 4)') 
            
        S=10**(-30.7667+26.49*np.log10(3328) - 7.0977*(np.log10(3328))**2 + 0.605*(np.log10(3328)**3) )
        
            
        
        file1=np.loadtxt(s,usecols=range(1,10))
    
        col1=file1[:,4]
        col2=file1[:,5]
        col3=file1[:,6]
        colav1=(col2+col3)/2
        
        yvalue1=[]
        xvalue1=[]
        yvalue4=colav1[np.where((col1>3335)&(col1<3336))]
        xvalue4=col1[np.where((col1>3335)&(col1<3336))]
        
        if n==3:
            yvalue4=colav1[np.where((col1>3335)&(col1<3336))]
            xvalue4=col1[np.where((col1>3335)&(col1<3336))]
        for i in range(len(col1)):
            if col1[i]>3335 and col1[i]<3336:
                if col1[i]>3335.5 and col1[i]<3335.62:
                    continue
                yvalue1.append(colav1[i])
                xvalue1.append(col1[i])
        if n==3:
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
        def func4(x, a, b,d,e,f):
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
        
        def monte_carlo_integration(x, y, y_err, n_simulations=1000):
            integrals = []
            for _ in range(n_simulations):
                y_simulated = y + np.random.normal(0, y_err, size=y.shape)
                integral = integrate.simpson(y_simulated, x)
                integrals.append(integral)
            integral_mean = np.mean(integrals)
            integral_std = np.std(integrals)
            return integral_mean, integral_std
        
        xvalue1=np.asarray(xvalue1)
        yvalue1=np.asarray(yvalue1) 
        p = np.polyfit(xvalue1, yvalue1, 6)
        error=[]
        for alpha in range(0,len(xvalue1)):
            
            error.append(np.std(yvalue1))
        
        #trial of a function to fit
        #if n==3:
            
            #popt, pcov = curve_fit(func4, xvalue1, yvalue1, maxfev=1000,sigma=error)
            #err=np.sqrt(abs(np.diag(pcov)))
        #if n==3:
            #y_model=func4(xvalue4,*popt)    
        
    
        y_model = polyN(xvalue4, p)
        
        
        
        
        
        
        xvalue3=col1[np.where((col1>3335.51)&(col1<3335.61))]
        yvalue3=colav1[np.where((col1>3335.51)&(col1<3335.61))]
        
        xvalue2=col1[np.where((col1>3335)&(col1<3336))]
        yvalue2=colav1[np.where((col1>3335)&(col1<3336))]
        
        if n==3:
            yvalue2=colav1[np.where((col1>3335.44)&(col1<3335.67))]
            xvalue2=col1[np.where((col1>3335.45)&(col1<3335.67))]
        
        
        convol1=convolve(yvalue4,Box1DKernel(1))
        
        flux=convol1-y_model
        
        
        convol1chopped=convol1[10:len(convol1)-10]
        xvalue2chop=xvalue2[10:len(xvalue2)-10]
        
        
        #test plots
        """
        plt.plot(xvalue4,yvalue4,label='with')
        plt.plot(xvalue1,yvalue1,label='wl') 
        plt.plot(xvalue4,y_model,label='model')
        """
        velocity=[]
        c=3*10**5
        for i in range(0,len(xvalue4)):
            v=((3335.4810-xvalue4[i])*c/xvalue4[i])
            velocity.append(v)
               
        velocity=np.asarray(velocity) 
        mask = ((xvalue4 > 3335) & (xvalue4 < 3335.5)) | ((xvalue4 > 3335.62) & (xvalue4 < 3336))
        velocitymask=velocity[mask]
        obser=yvalue4[mask]
        
        
        model=y_model[mask]
        final=flux[mask]
        
        modelx=xvalue4[mask]
        xvaluemask=xvalue4[mask]
        mask2 = ((xvalue4 > 3335.5) & (xvalue4 < 3335.62))
        yint=flux[mask2]
        
        
        NCH2=integrate.simpson(yint,dx=velocity[0]-velocity[1])
        baseline=[]
        add=0
        for beta in range(0,len(flux)):
            add=flux[beta]**2+add
        
        if n==3:
#trial of basline error      
                #for zeta in range(0,len(final)):
                #baseline_unc=np.sqrt(model[zeta]**2*((err[0]/popt[0])**2)+(err[1]/popt[1])**2+(err[2]/popt[2])**2+(err[3]/popt[3])**2+(err[4]/popt[4])**2)
                #baseline.append(baseline_unc)
        
            sigma3=np.sqrt(np.std(final)**2+add/np.sqrt(len(flux)))
        else:
            sigma3=np.std(final)
        
        integral3, integral3_unc = monte_carlo_integration(velocity, flux, sigma3)
        
        I.append(NCH2)
        errorval.append(integral3_unc)
        
        plt.title("3335 MHz spectral line",fontsize=20)

        colour=['green','blue','red']
        plt.plot(velocity,flux,color=colour[n-2],linestyle='-',label=po)
        if n==3:
            plt.axhline(y=0,color='k',label='Zero Flux line')
        plt.xlabel("Frequency (MHz)",fontsize=20)
        plt.ylabel("Flux (Jansky)",fontsize=20)
        plt.legend(fontsize=20)
        plt.xlim(-40,30)
        
with open('int_intensityCH.txt', 'w') as wf:
     writer=csv.writer(wf,delimiter=' ')
     
     Column1=I
     Column2=errorval
     writer.writerows(zip(Column1,Column2))
#Gaussian fit
"""
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

    
"""            
#histogram       
"""



yvaluesubhist=yvalue3-func2(xvalue3, *popt)

suma=0

for beta in range(0,len(yvaluesubhist)):
    
    suma=suma+(yvaluesubhist[beta])**2
    

alpha=integrate.simpson(yvaluesubhist,dx=velocity[0]-velocity[1])
print(alpha)


plt.xlabel("Frequency (MHz)",fontsize=20)
plt.ylabel("Flux (Jansky)",fontsize=20)
plt.legend(fontsize=20)

hist,binedge=np.histogram(yvaluesubhist,bins=10,density=False)
weights=np.ones_like(yvaluesubhist)/len(yvaluesubhist)

plt.hist(yvaluesubhist, weights=weights)

  
plt.title("3335 MHz spectral line",fontsize=20)
plt.savefig('CHmet4',dpi=300)
"""