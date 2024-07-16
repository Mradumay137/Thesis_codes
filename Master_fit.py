#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 00:00:36 2024

@author: maverick
"""

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
import csv
from scipy.stats import norm
import statistics
from collections import Counter
from scipy.fft import fft, ifft

rfi1667=[]
rfi1665=[]
NOHvalues1=[]
NOHvalues2=[]
stdev=[]
sigmaarr1667=[]
sigmaarr1665=[]
meanarr=[]
badbase=[]
position=[]
for u in range(1,3):
    baseline=[]
    itera=[]
    stdmn=[]
    stdnp=[]
    for j in range(1,52):
        if j==50:
            continue
        print("Processing file {0}.. Please wait".format(j))
        itera.append(j)
        s=f'S{j}.dat'
        if j==41:
            s='E1.dat'
        if j==42:
            s='E2.dat'
        if j==43:
            s='E3.dat'
        if j==44:
            s='W1.dat'
        if j==45:
            s='W2.dat'
        if j==46:
            s='W3.dat'    
        if u==1:
            C1=2.2
            fmin=1667
            fmax=1668
            flinemin=1667.38
            flinemax=1667.44
            fres=1667.3590
            vmin=-25
            vmax=10
            Tex=4
        if u==2:
            C1=4
            fmin=1665.2
            fmax=1665.6
            flinemin=1665.42
            flinemax=1665.46
            fres=1665.4019
            vmin=-18
            vmax=0
            Tex=3.72
        if j==47 or j==48 or j==49 or j==51:
            s=f'S{j-41}OH.dat'
               
        
        file=np.loadtxt(s,usecols=range(1,10),skiprows=2)
        
        col1=file[:,4]
        col2=file[:,5]
        col3=file[:,6]
        colav1=(col2+col3)/2
         
        if j==47 or j==48 or j==49 or j==51:
            col1=file[:,4]
            col2=file[:,5]
            col3=file[:,6]
            colav1=(col2+col3)/2 
        
        
        col4=[]
        col5=[]
        
        xwoline=[]
        ywoline=[]
        def func0(x,m,inter):
            return m*x+inter
        
        def func1(x, a, b, c,d,e):
            return a*np.sin(x)**2 + b*np.cos(x)**2+c*np.sin(x)+d*np.cos(x)+e
        
        
        def func2(x, a, b, c,d,e):
            return a*x**4+b*x**3+c*x**2+d*x+e
        
        def func3(x, a, b,c,d):
            return a*x**3+b*x**2+c*x+d
        
        xwline=[]
        ywline=[]
        for i in range(len(col1)):
            if col1[i]>fmin and col1[i]<fmax:
                xwline.append(col1[i])
                ywline.append(colav1[i])
            
            if col1[i]>fmin and col1[i]<fmax:
                if col1[i]>flinemin and col1[i]<flinemax:
                    continue
                ywoline.append(colav1[i])
                xwoline.append(col1[i])
       
        popt2, pcov2 = curve_fit(func2, xwoline, ywoline)
        
        ywline=np.asarray(ywline)
        xwline=np.asarray(xwline)
        
        #RFI removal
        
        """
        
        gamma=100
            
        x=np.linspace(fmin,fmax,gamma)
        for beta in range(0,gamma-2):
            
            low=x[beta]
            high=x[beta+1]
            
            mask=(xwline > low) & (xwline < high)
            ybin= ywline[mask]
            xbin = xwline[mask]
            dx=xwline[1]-xwline[0]
            
            dy=np.gradient(ybin, dx)
            
            if np.std(dy)>20:
                exclusion_mask = ~np.isin(ywline, ybin)
                ywline = ywline[exclusion_mask]
                xwline = xwline[exclusion_mask]
               
        """
        
       
             
        def polyN(x, p):
        	y = 0
        	i = len(p) - 1
        	p = np.flip(p)
        	while (i >= 0):
        		#print(i)
        		y += p[i]*x**float(i)
        		i -= 1
        	return y
        
        
        
        ywoline=np.asarray(ywoline)
        xwoline=np.asarray(xwoline)
        
        
       
        p = np.polyfit(xwoline, ywoline, 6)
        
        
        
        
        y_model = polyN(xwline, p)
        
        
        model_rot=polyN(xwoline, p)
        
        
            
            
            
        
        ywline=np.asarray(ywline)
        xwline=np.asarray(xwline)
        
        
        
        #ysub=ywline-func2(xwline, *popt2)
        
        ysub2=ywline-y_model
        
        
        
        """
        mask=((xwline > flinemin-0.06) & (xwline < flinemax+0.06)) 
        ymask=ysub[mask]
        ysub2=ysub-abs(np.mean(ymask))
        """
        
        
    
        
        #popt3,pcov3=curve_fit(func0,xwline,ysub)
        
        #if u==2:
            #ysub=ywline-y_model
        
        #ysub2=ysub-func0(xwline,*popt3)
        
        
             
        
        
        #check efficacy of baseline fit
        
        """
        mean=np.mean(ysub2[mask])
        if mean>0.0005:
            badbase.append(mean)
            position.append(j)
        baseline.append(mean)
        """
        
            
        def func(x,a,b,c):
            return a*np.exp(-(x-b)**2/(2*c**2))
        
        """
        poptg,pcovg=curve_fit(func,xwline,ysub)
        
        
        mask=((xwline > flinemin-0.08) & (xwline < flinemin)) | ((xwline > flinemax) & (xwline < flinemax+0.08))
        ywoline2=ysub[mask]
        xwoline2=xwline[mask]
        p3= np.polyfit(xwoline2, ywoline2, 10)
        y_model_new=polyN(xwline,p3)
        ysub=ysub-y_model_new
        """    
        for i in range(0,len(xwline)):
            if xwline[i]>fmin and xwline[i]<fmax:
                col4.append(xwline[i])
                col5.append(ysub2[i])
        
        velocity=[]
        cs=3*10**5
        for i in range(0,len(col4)):
            v=((fres-col4[i])*cs/col4[i])
            velocity.append(v)
        
        col5smooth=convolve(col5,Box1DKernel(1))
        velocity=np.asarray(velocity)
        
        if j>=1 and j<=3:
            vmin=-10
            vmax=-4
        
        if j>=4 and j<=5:
            vmin=-15
            vmax=-6
        if j>5 and j<9:
            vmin=-11
            vmax=-4
        if j==9 or j==10:
            vmin=-13
            vmax=-3
        if j>10 and j<15:
            vmin=-10
            vmax=-4
        if j==15:
            vmin=-12 
            vmax=-6
            
        if j>15 and j<20:
            vmin=-10 
            vmax=-4
        if j>=20 and j<=30:
            vmin=-8
            vmax=-1
            
        if j==31 or j==32:
            vmin=-10
            vmax=0
        
        if j==33:
            vmin=-13
            vmax=-3
        if j==40:
            vmin=-16
            vmax=-0
        if j==34:
            vmin=-9
            vmax=-5
        if j>34 and j<39:
            vmin=-8
            vmax=0
        if j==39:
            vmin=-6
            vmax=0
        if j==41:
            vmin=-10
            vmax=-4
        if j==42:
            vmin=-8
            vmax=-2
            
        if j==43:
            vmin=-7.5
            vmax=-3
        if j>43 and j<45:
            vmin=-12
            vmax=-2
        if j==45 or j==46:
            vmin=-14
            vmax=-8 
        if j>46 and j<52:
            vmin=-10
            vmax=-3
        
        mask2=((velocity > vmin-5) & (velocity < vmax+5)) 
        velocitysmall=velocity[mask2]
        col5smoothsmall=col5smooth[mask2]
        
        psmall=np.polyfit(velocitysmall,col5smoothsmall, 8)
        
        modelmask=polyN(velocitysmall, psmall)
        
        
        
        plt.figure(figsize=(5,5))
        
        plt.plot(velocitysmall,modelmask,label=j)
        
        
        
        plt.plot(velocitysmall,col5smoothsmall)
        
        plt.legend()
        
        
        
        
        def first_derivative(array, dx):
            return np.gradient(array, dx)

        def second_derivative(array, dx):
            first_deriv = np.gradient(array, dx)
            return np.gradient(first_deriv, dx)

        dx=xwline[1]-xwline[0]
        fderi = first_derivative(modelmask, dx)
        sderi = second_derivative(modelmask, dx)
        mid_index = len(velocitysmall) // 2

        R = abs(((1 + fderi[mid_index]**2)**1.5) / sderi[mid_index])
        
        print(R)
        
        bins=[velocity[0],velocity[1],velocity[2]]
        for eta in range(3,len(velocity)):
            bin1=(velocity[eta]+velocity[eta-1]+velocity[eta-2]+velocity[eta-3])/4
            bins.append(bin1)
        bins=np.asarray(bins)
        popt, pcov = curve_fit(func, bins,col5smooth,p0=(0.04,-7,0.2),method='lm',maxfev=10000)
        if j==39:
            popt, pcov = curve_fit(func, bins,col5smooth,p0=(0.04,-3,0.2),method='lm',maxfev=10000)
        if j==37:
            popt, pcov = curve_fit(func, bins,col5smooth,p0=(0.04,-4,0.2),method='lm',maxfev=10000)
        if j==45:
            popt, pcov = curve_fit(func, bins,col5smooth,p0=(0.04,-10,0.2),method='lm',maxfev=10000)
        if j==46:
            popt, pcov = curve_fit(func,bins,col5smooth,p0=(0.04,-12,0.2),method='lm',maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
        xerr=[]
        for alpha in range(0,len(xwline)):
            
            xe=cs*0/xwline[i]
            xerr.append(xe)
            
    
        a=popt[0]
        b=popt[1]
        c=popt[2]
        """
        d=popt[3]
        e=popt[4]
        f=popt[5]
    
        """
        """
        Nstd=[]
        for i in range(0,len(xwline)):
            sigma = np.exp(-(bins[i] - b)**2 / (2 * c**2)) * ((perr[0]**2) - a * ((bins[i] - b) / c)**2 * (perr[1]**2) - 2 * a * (perr[2]**2) / c**3 + a * (bins[i]- b) / c**2 * (xerr[i])**2)
            sig=math.sqrt(abs(sigma))
            Nsd=sig*C1*10**14*4
            Nstd.append(Nsd)
        Nsig=np.std(Nstd)
        """
        
            
        
        """
        
        """
        bgsd=np.std(ywoline)
        linearr=[]
        varr=[]
        if u==1:
            for gamma in range(0,len(col5smooth)):
        
                if xwline[gamma]>flinemin and xwline[gamma]<flinemax:
                    linearr.append(col5smooth[gamma])
                    varr.append(velocity[gamma])
            maxi=np.argmax(linearr)        
            vcentral=varr[maxi]
            
           
            
        if u==2:
            for gamma in range(0,len(col5smooth)):
                if xwline[gamma]>flinemin and xwline[gamma]<flinemax:
                    linearr.append(col5smooth[gamma])
                    varr.append(velocity[gamma])
            maxi=np.argmax(linearr)        
            vcentral=varr[maxi]
            
           
        
        
    
        col6=col5smooth[np.where((velocity>vmin)&(velocity<vmax))]
        velocitychop=bins[np.where((velocity>vmin)&(velocity<vmax))]
        
        
        
        #popt, pcov = curve_fit(func, velocitychop,col6,p0=(0.04,-(abs(vmin)-abs(vmax))/2,0.2),method='lm',maxfev=10000)
        
        """
        plt.figure(figsize=(5,5))
        plt.plot(velocity,col5smooth,label=j)
        plt.xlim(-25,25)
        
        plt.legend()
        """
        """
        plt.figure(figsize=(5,5))
        plt.plot(velocitychop,col6,label=j)
        plt.figure(figsize=(5,5))
        plt.plot(velocity,col5smooth)
        plt.axvline(x=vmin)
        plt.axvline(x=vmax)
        plt.legend()
        """
        N=integrate.simpson(col6,dx=velocitychop[0]-velocitychop[1])
        
        #plt.figure(figsize=(5,5))
        #plt.plot(velocitychop,col6,label=j)
        #plt.legend()
        NOH=abs(N)*(C1*10**(14))*(Tex/(Tex-3))*0.87
        
        if u==1:
            if R>0.051:
                NOH=0
        
        if u==2:
        
            if R > 0.1:
                NOH=0
        
        mask=((velocity > 10) & (velocity < 16)) 
        sigcol=col5smooth[mask]
        sigma3=np.std(sigcol)
        #sigma3=integrate.simpson(sigcol,dx=velocitychop[0]-velocitychop[1])
        Nsig=(velocitychop[0]- velocitychop[1])*np.sqrt(len(col6))*sigma3*C1*10**14*Tex/(Tex-3)
        
        if j==1:
            Nsigref=Nsig
                
        """
        plt.figure(figsize=(5,5))
        plt.plot(velocity,col5smooth,'k',label=j)
        plt.legend()
        plt.xlim(-25,25)
        """
        """
        plt.axhline(y=0,color='black',linestyle='dashdot')
        if u==1:
            plt.legend()
            plt.savefig('S{0} 1667'.format(j))
        if u==2:
            plt.legend()
            plt.savefig('S{0} 1665'.format(j))     
        plt.close()
        """
        """
        if NOH<3*10**13:
            Nsig=0
        """
        if u==1:
            sigmaarr1667.append(Nsig)
        if u==2:
            sigmaarr1665.append(Nsig)
        
        
        if u==1:
            NOHvalues1.append(NOH)
        if u==2:
            NOHvalues2.append(NOH)
        
        if u==1:
            rfi1667.append(Nsig)
        if u==2:
            rfi1665.append(Nsig)
            
    
        #plt.figure(figsize=(5,5))
        """
        hist,binedge=np.histogram(col5smooth,bins=10,density=False)
        hist=hist/sum(hist)
        
        n2=len(hist)
        x2=np.zeros((n2),dtype=float)
        for kappa in range(n2):
            x2[kappa]=(binedge[kappa+1]+binedge[kappa])/2
        y=hist
        
        xmin,xmax=plt.xlim()
        
       
        mean = sum(x2*y)/sum(y) 
        meanarr.append(mean)                 
        sigma = sum(y*(x2-mean)**2)/sum(y) 
        
        print(sigma)
        
        popt5,pcov5=curve_fit(func,x2,y, p0=(max(y),mean, 0.5),maxfev=5000)
        
        
        #plt.figure(figsize=(10,7))
        x1=np.linspace(-0.05,0.05,1000)
        
        
        weights=np.ones_like(col5smooth)/len(col5smooth)
        
        """
        
        #plt.hist(col5smooth, weights=weights,label='Data')
        
        
        #plt.xlim(-0.05,0.05)
        
         
        #plt.plot(x1, func(x1,*popt5), 'k', linewidth=2, label='Gaussian fit')
        
        #plt.xlabel('Flux density values (Jy)',fontsize=18)
        
        #plt.ylabel('Normalized frequency',fontsize=18)
        
        #plt.legend(fontsize=14)
        #plt.savefig('S{0} histogram'.format(j),dpi=200)
        
        #plt.show()
        
        #plt.close()
        
        """
        
        NOHvalues.append(NOH)
        
        
        
        if j==1:
            
        if j==2:
            
    
        plt.legend()
        plt.axhline(y=0,color='black',linestyle='dashdot')
        plt.xlabel('Velocity (Km/s)',fontsize=14)
        plt.ylabel('Flux (Jansky)',fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.savefig("Fit",dpi=200,bbox_inches='tight')
        plt.legend()
        """
        #plt.figure(figsize=(5,5))
        #plt.plot(bins,col5smooth,'r',label=s)
        #plt.xlim(-25,25)
        
        
        
        
        #plt.plot(velocity, func(velocity,*popt),'b',label='Fit')
        
        #plt.axvline(x=vmin)
        #plt.axvline(x=vmax)
        
        
        #plt.axhline(y=0)
        """
        if u==1:
            plt.legend()
            plt.savefig('S{0} 1667'.format(j))
        if u==2:
            plt.legend()
            plt.savefig('S{0} 1665'.format(j))     
        plt.close()
        """
    if u==1:        
        NOHvalues1[5] = (NOHvalues1[5]+NOHvalues1[46])/2
        NOHvalues1[6] = (NOHvalues1[6]+NOHvalues1[47])/2
        NOHvalues1[7] = (NOHvalues1[7]+NOHvalues1[48])/2
        NOHvalues1[9] = (NOHvalues1[9]+NOHvalues1[49])/2
        del NOHvalues1[49]
        del NOHvalues1[48]
        del NOHvalues1[47]
        del NOHvalues1[46]
    if u==2:
        NOHvalues2[5] = (NOHvalues2[5]+NOHvalues2[46])/2
        NOHvalues2[6] = (NOHvalues2[6]+NOHvalues2[47])/2
        NOHvalues2[7] = (NOHvalues2[7]+NOHvalues2[48])/2
        NOHvalues2[9] = (NOHvalues2[9]+NOHvalues2[49])/2
        del NOHvalues2[49]
        del NOHvalues2[48]
        del NOHvalues2[47]
        del NOHvalues2[46]
       
for omicron in range(1,len(NOHvalues1)):
    
    if abs(NOHvalues1[omicron])> abs(1.5*NOHvalues2[omicron]) and abs(NOHvalues1[omicron])>0.5*10**14 and abs(NOHvalues2[omicron])>0.5*10**14 : 
        #NOHvalues1[omicron]=NOHvalues2[omicron]
        print(omicron,NOHvalues1[omicron]/10**14,NOHvalues2[omicron]/10**14)
        
     
    if abs(NOHvalues2[omicron])> abs(1.5*NOHvalues1[omicron]) and abs(NOHvalues2[omicron])>0.5*10**14 and abs(NOHvalues1[omicron])>0.5*10**14:
        #NOHvalues2[omicron]=NOHvalues1[omicron]
        print(omicron,NOHvalues1[omicron]/10**14,NOHvalues2[omicron]/10**14)
  
with open('Output.txt', 'w') as wf:
     writer=csv.writer(wf,delimiter=' ')
     Column1=NOHvalues1
     Column2=NOHvalues2
     Column3=sigmaarr1667
     Column4=sigmaarr1665
     writer.writerows(zip(Column1,Column2,Column3,Column4))

       
   