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
import cblind as cb
from scipy.integrate import quad
from scipy import integrate
from scipy.fft import fft, fftfreq
import csv
    
plt.figure(figsize=(10,10))
sub=0
meanarr=[]
integral_uncertainties = []
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
I=[]
error=[]
for omicron in range(6,11):
    if omicron==9:
        continue
    
    plt.figure(figsize=(10,10))
    for n in range(1,8,2):
       
        if n==1:
            #s=str('2a_scal_new.dat')
            s=str(f'S{omicron}met1.dat')
            #s=str('2a_scal_new.dat')
        if n==2:
            s=str('2b_tcal.dat')  
        if n==3:
            #s=str('method2wns.dat')
            s=str(f'S{omicron}met2.dat')
        if n==4:
            s=str('2c_tcal.dat')
        if n==5:
            #s=str('method3.dat')
            s=str(f'S{omicron}met3.dat')
        if n==6: 
            s=str('2d_tcal.dat')
        if n==7:
            
            #s=str('S6OH.dat')
            s=str(f'S{omicron}OH.dat')
            #s=str('method5.dat')
        if n==1:
            po=str('ON-OFF/OFF (Method 1)')
        if n==2:
           plt.title("Method 1 (Kelvin)",fontsize=14) 
        if n==3:
           po=str('ON-OFF with noise cal (Method 2)')  
        if n==4:
            plt.title("Method 2 (Kelvin)",fontsize=14)
        if n==5:
           po=str('ON with noise cal (Method 3)')
        if n==6:
            plt.title("Method 3 (Kelvin)",fontsize=14)
        if n==7:
            po=str('ON with 1934 (Method 4)')
               
            
        
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
            
            if col1[i]>1665.3 and col1[i]<1665.6:
                yvalueoriginal2.append(colav1[i]) 
                xvalueoriginal2.append(col1[i])
            
            if n==7:
                
                if col1[i]>1667.2 and col1[i]<1667.6:
                    yvalueoriginal3.append(colav1[i]) 
                    xvalueoriginal3.append(col1[i])
                    
            elif n==5:
                if col1[i]>1667.2 and col1[i]<1667.6:
                    yvalueoriginal3.append(colav1[i]) 
                    xvalueoriginal3.append(col1[i])
                    
            else:
                if col1[i]>1667.1 and col1[i]<1667.7:
                    yvalueoriginal3.append(colav1[i]) 
                    xvalueoriginal3.append(col1[i])
            
            if col1[i]>1612.1 and col1[i]<1612.5:
                yvalue1.append(colav1[i])
                xvalue1.append(col1[i])
            if col1[i]>1665.3 and col1[i]<1665.6:
                if col1[i]>1665.42 and col1[i]<1665.46:
                    continue
                yvalue2.append(colav1[i])
                xvalue2.append(col1[i])
                
            
                
            if n==7 :
                
                if col1[i]>1667.2 and col1[i]<1667.6:
                    if col1[i]>1667.38 and col1[i]<1667.45:
                        continue
                    if col1[i]>1667.73 and col1[i]<1667.80:
                        continue
                    yvalue3.append(colav1[i])
                    xvalue3.append(col1[i])
            elif n==5:
                if col1[i]>1667.2 and col1[i]<1667.6:
                    if col1[i]>1667.38 and col1[i]<1667.45:
                        continue
                    if col1[i]>1667.73 and col1[i]<1667.80:
                        continue
                    yvalue3.append(colav1[i])
                    xvalue3.append(col1[i])
            
            else:
                if col1[i]>1667.1 and col1[i]<1667.7:
                    
                    if col1[i]>1667.38 and col1[i]<1667.45:
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
        delarr=[]
        
        
        c=3*10**5
        
        
        
            
        
        
        def monte_carlo_integration(x, y, y_err, n_simulations=1000):
            integrals = []
            for _ in range(n_simulations):
                y_simulated = y + np.random.normal(0, y_err, size=y.shape)
                integral = integrate.simpson(y_simulated, x)
                integrals.append(integral)
            integral_mean = np.mean(integrals)
            integral_std = np.std(integrals)
            return integral_mean, integral_std
            
        def func2(x, a, b, c,d):
            return a*x**3+b*x**2+c*x+d
        
        def func(x, a, b):
            return a*x+b
            """
            return a*x**4+b*x**3+c*x**2+d*x**1+e
            """
        def func5(x,al,bl,cl):
            return al*np.exp(-(x-bl)**2/(2*cl**2))
        
       
        
        def polyN(x, p):
        	y = 0
        	i = len(p) - 1
        	p = np.flip(p)
        	while (i >= 0):
        		#print(i)
        		y += p[i]*x**float(i)
        		i -= 1
        	return y
        xvalue3=np.asarray(xvalue3)
        yvalue3=np.asarray(yvalue3) 
        xvalueoriginal3=np.asarray(xvalueoriginal3)
        yvalueoriginal3=np.asarray(yvalueoriginal3)
        #p, p_err = np.polyfit(xvalue3, yvalue3, 6)
        
        #RFI flagging
        prep=np.polyfit(xvalue3, yvalue3, 2)
        y_model_pre = polyN(xvalueoriginal3,prep)
        y_model_pre_wl=polyN(xvalue3,prep)
        for beta in range(0,len(yvalue3)):
            if abs(yvalue3[beta]-y_model_pre_wl[beta])>1.5*np.std(yvalue3-y_model_pre_wl) and xvalue3[beta]>1667.43:
                delarr.append(beta)
            elif abs(yvalue3[beta]-y_model_pre_wl[beta])>1.5*np.std(yvalue3-y_model_pre_wl) and xvalue3[beta]<1667.38:  
                delarr.append(beta)
        xvalue3=np.asarray(xvalue3)
        yvalue3=np.asarray(yvalue3)   
        yvalue3=np.delete(yvalue3,delarr)       
        xvalue3=np.delete(xvalue3,delarr)
        
        delarr2=[]
        for beta in range(0,len(yvalueoriginal3)):
            if abs(yvalueoriginal3[beta]-y_model_pre[beta])>1.5*np.std(yvalueoriginal3-y_model_pre) and xvalueoriginal3[beta]>1667.43:
                delarr2.append(beta)
            elif abs(yvalueoriginal3[beta]-y_model_pre[beta])>1.5*np.std(yvalueoriginal3-y_model_pre) and xvalueoriginal3[beta]<1667.38:  
                delarr2.append(beta)
         
        yvalueoriginal3=np.delete(yvalueoriginal3,delarr2)       
        xvalueoriginal3=np.delete(xvalueoriginal3,delarr2)
        #RFI flagging
        
        convol1=convolve(yvalue1,Box1DKernel(1))
        convol1chopped=convol1[10:len(convol1)-10]
        
        convol2=convolve(yvalueoriginal2,Box1DKernel(1))
        convol2chopped=convol2[10:len(convol2)-10]
        
        convol3=convolve(yvalueoriginal3,Box1DKernel(1))
        convol3chopped=yvalueoriginal3[10:len(convol3)-10]
        
        convol4=convolve(yvalue4,Box1DKernel(1))
        convol4chopped=convol4[10:len(convol4)-10]
        
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
        
        xvalue1=np.asarray(xvalue1)
        xvalue2=np.asarray(xvalue2)
        xvalue3=np.asarray(xvalue3)
        yvalueoriginal3=np.asarray(yvalueoriginal3)
        xvalue4=np.asarray(xvalue4)
        xvalueoriginal2=np.asarray(xvalueoriginal2)
        xvalueoriginal3=np.asarray(xvalueoriginal3)
        xvalueoriginalchop2=xvalueoriginal2[10:len(xvalueoriginal2)-10] 
       
        """
        if n==5 or n==6 :
            popt2, pcov2 = curve_fit(func2, xvalue2, yvalue2)
            
            yvaluesub2=convol2chopped-func2(xvalueoriginalchop2, *popt2)
        else:
            popt2, pcov2 = curve_fit(func2, xvalue2, yvalue2)
            yvaluesub2=convol2chopped-func2(xvalueoriginalchop2, *popt2)
        
        if n==5 or n==6 :
            popt1, pcov1 = curve_fit(func2, xvalue1, yvalue1)
            xvalue1=xvalue1[10:len(xvalue1)-10]
            yvaluesub1=convol1chopped-func2(xvalue1, *popt1)
        else:
            popt1, pcov1 = curve_fit(func2, xvalue1, yvalue1)
            xvalue1=xvalue1[10:len(xvalue1)-10]
            yvaluesub1=convol1chopped-func2(xvalue1, *popt1)
        
        """
        #convert to 1665 line (next 4 lines)
        """
        xvalue3=xvalue2
        yvalue3=yvalue2
        xvalueoriginal3=xvalueoriginal2
        convol3chopped=convol2chopped
        """
        xvalueoriginalchop3=xvalueoriginal3[10:len(xvalueoriginal3)-10]
        
        p= np.polyfit(xvalue3, yvalue3, 6)
        
        
            
        
        
    
        y_model = polyN(xvalueoriginalchop3,p)
        
        
        
        
        
        
        popt7, pcov7 = curve_fit(func, xvalue3, yvalue3)
        yvaluesub3=convol3chopped-y_model
        
        
        
        if n==1 and omicron==6:
            yvaluesub3=yvaluesub3*0.986 #because the conversion factor is 35 not 35.5
           
        
        popt4, pcov4 = curve_fit(func2, xvalue4, yvalue4)
        xvalue4=xvalue4[10:len(xvalue4)-10]
        yvaluesub4=convol4chopped-func2(xvalue4, *popt4)
        
        def func4(x,a,b,c,d,e,f):
            return a*np.exp(-(x-b)**2/(2*c**2))+d*np.exp(-(x-e)**2/(2*f**2))
        
    
        velocitychop1=velocity1[10:len(velocity1)-10]
        velocitychop1=np.asarray(velocitychop1)
        velocitychop2=velocity2[10:len(velocity2)-10]
        velocitychop2=np.asarray(velocitychop2)
        velocitychop3=velocity3[10:len(velocity3)-10]
        velocitychop3=np.asarray(velocitychop3)
        velocitychop4=velocity4[10:len(velocity4)-10]
        velocitychop4=np.asarray(velocitychop4)
        
        #velocitychop3=velocitychop2 # convert to 1665
        
        #plt.plot(velocitychop3,yvaluesub3,label=po)
        
    
        #plt.axhline(y=0,color='black',linestyle='dashdot',label='Flux=0 line')
        #plt.legend()
        
        
        
        
        # following line is for 1665 
        #mask = ((xvalueoriginalchop3 > 1665.2) & (xvalueoriginalchop3 < 1665.42)) | ((xvalueoriginalchop3 > 1665.46) & (xvalueoriginalchop3 < 1665.7))
        
        #1667 Line
        
        mask = ((xvalueoriginalchop3 > 1667.1) & (xvalueoriginalchop3 < 1667.38)) | ((xvalueoriginalchop3 > 1667.45) & (xvalueoriginalchop3 < 1667.7))
        if n==1 or n==3:
            mask = ((xvalueoriginalchop3 > 1667.1) & (xvalueoriginalchop3 < 1667.38)) | ((xvalueoriginalchop3 > 1667.45) & (xvalueoriginalchop3 < 1667.7))
        yvaluesubsampl = yvaluesub3[mask]
        
        
        #testing frequency range
        """
        xvaluesampl=velocitychop3[mask]
        xvaluesamplfreq=xvalueoriginalchop3[mask]
        model=y_model[mask]
        print(np.std(yvaluesubsampl))
        if n==1 or n==3:
            plt.plot(xvaluesamplfreq,yvaluesubsampl)
            plt.plot(xvalue3,yvalue3-np.mean(yvalue3))
            plt.axhline(y=0,color='black',linestyle='dashdot',label='Flux=0 line')
            plt.plot(xvaluesamplfreq,model)
            plt.axhline(y=np.mean(yvalue3),color='black',linestyle='dashdot',label='Flux=0 line')
        
        """
        #Gaussian fit
        """
        popt, pcov = curve_fit(func4, velocitychop3,yvaluesub3,p0=(0.04,-7,0.2,0.05,-8,0.2),method='lm',maxfev=5000)
        
        a=popt[0]
        b=popt[1]
        c=popt[2]
        d=popt[3]
        e=popt[4]
        f=popt[5]
        
        #NOH=quad(func4,velocitychop3[-1],velocitychop3[0],args=(a,b,c,d,e,f))
        """
        #1667
        mask2 = ((xvalueoriginalchop3 > 1667.37) & (xvalueoriginalchop3 < 1667.45))
        #1665  
        #mask2 = ((xvalueoriginalchop3 > 1665.42) & (xvalueoriginalchop3 < 1665.46))
        yvalueline=yvaluesub3[mask2]
        xvalueline=xvalueoriginalchop3[mask2]
        NOH2=integrate.simpson(yvalueline,dx=velocitychop3[0]-velocitychop3[1])
        
        su=np.sum(yvaluesub3**2)
        
        
        print(NOH2)
        I.append(NOH2)
        sigma3 = np.std(yvaluesubsampl)
        
        if n==1:
            sigma3=np.sqrt(sigma3**2+su/(4*35**2))
        integral3, integral3_unc = monte_carlo_integration(velocitychop3, yvaluesub3, sigma3)
        
        print(integral3_unc)
        error.append(integral3_unc)
        #baseline_fit_error
        """
        p,cov= np.polyfit(xvalue3, yvalue3, 2,cov=True,full=False)
        if n==1 or n==3:
            p,cov= np.polyfit(xvalue3, yvalue3, 1,cov=True,full=False)
        
        err=np.sqrt(abs(np.diag(cov)))
        print(err)
        print(p)
        baseline=[]
        for zeta in range(0,len(yvaluesubsampl)):
            baseline_unc=err[0]**2* (model)**2+err[1]**2
            baseline.append(baseline_unc)
        print(np.std(baseline))
        """
        """
        color,linestyle=cb.Colorplots().cblind(5)
        sigma3 = np.std(yvaluesub3)
        integral3, integral3_unc = monte_carlo_integration(velocitychop3, yvaluesub3, sigma3)
        plt.plot(velocitychop3,yvaluesub3,label=po)
        if n==7:
            plt.axhline(y=0,color='black',linestyle='dashdot',label='Flux=0 line')
        plt.legend()
        """
        #standard error
        #print((velocitychop3[0]-  velocitychop3[1])*np.sqrt(len(yvaluesub3))*sigma3)
        
    
        
        
        
    
        
        
        #Baseline_fits
        """
        if n==5 or n==7:
            convol3chopped = convol3chopped - 33.7
            y_model= y_model - 33.7
        if n==3:
            convol3chopped = convol3chopped - 0.1
            y_model= y_model - 0.1
            
       
        if n==7:
            ax2 = plt.subplot(gs[1])
            ax2.plot(velocitychop3,convol3chopped,label=po)
            ax2.plot(velocitychop3,y_model)
            plt.xlabel("Velocity (km/s)",fontsize=20)
            plt.ylabel("Flux (Jansky)",fontsize=2)
            plt.xlim(-40,20)
        else:
            ax1 = plt.subplot(gs[0])
            ax1.plot(velocitychop3,convol3chopped,label=po)
            plt.plot(velocitychop3,y_model)
            plt.xlim(-40,20)
        plt.legend(fontsize=14)
        
        plt.xlabel("Velocity (km/s)",fontsize=20)
        plt.ylabel("Flux (Jansky)",fontsize=20)
        plt.title("1667",fontsize=20)
        plt.savefig('Baseline_fits', dpi=300) 
        
        plt.xlim(-40,20)
        """
        #spectral line plots with offset
    
       
        
        
        
        
       
        plt.plot(velocitychop3,yvaluesub3+sub,label=po)
        sub=sub+0.0
        plt.legend(fontsize=14)
        plt.xlabel("Velocity (km/s)",fontsize=18)
        plt.ylabel("Flux (Jansky)",fontsize=18)
        
       
         
         
            
        
        if n==5:
            plt.axhline(y=0,color='black',linestyle='dashdot',label='Flux=0 line')
       
        plt.xlim(-40,30)
        
       
       
        
        
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.savefig('Juxtaposed1934ns', dpi=300)
        
        #histogram code
        
        """
        hist,binedge=np.histogram(yvaluesubsampl,bins=21,density=False)
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
        
        
        
        #popt5,pcov5=curve_fit(func5,x2,y, p0=(0.22,0, 0.01))
        
        x1=np.linspace(-0.05,0.05,1000)
        
        
        weights=np.ones_like(yvaluesubsampl)/len(yvaluesubsampl)
        
       
        
        plt.hist(yvaluesubsampl, bins=11, weights=weights,label=po)
        
        
        plt.xlabel('Flux density values (Jy)',fontsize=18)
        
        plt.ylabel('Normalized number of values',fontsize=18)
        
        plt.legend(fontsize=14)
        
        plt.savefig('All_methods_histogram',dpi=300)
        
        plt.plot(x1, func5(x1,*popt5), 'k', linewidth=2, label='Gaussian fit')
        plt.show()
        
        """
        
        #FFT
        """
        ft=np.fft.fft(yvalue3)
        xf = fftfreq(len(yvalue3), 1 / velocitychop3[-1])
        
        plt.plot(xf,abs(ft))
        plt.ylim(0,5)
        plt.xlim(5,20)
        """
    
# Open the file in write mode
with open('int_intensity_corr_1667.txt', 'w') as wf:
     writer=csv.writer(wf,delimiter=' ')
     
     Column1=I
     Column2=error
     writer.writerows(zip(Column1,Column2))
"""
with open('int_intensity.txt', 'w') as file:
    for item in I:
        file.write(f"{item}\n")

    """        
    #other lines
    
"""
    plt.title("1665",fontsize=20)
    plt.plot(velocitychop2,yvaluesub2,label=p)
    plt.legend(fontsize=20)
    plt.xlim(-40,30)
    plt.ylim(-0.035,0.1)
    plt.xlabel("Velocity (km/s)",fontsize=22)
    plt.ylabel("Flux (Jansky)",fontsize=22)
    if n==5:
        plt.axhline(y=0,color='black',linestyle='dashdot',label='Flux=0 line')
    plt.subplot(7, 1, n)
    
    plt.plot(velocitychop1,yvaluesub1,label='1612')
    
    plt.subplot(2, 1, 2)
    plt.plot(velocitychop2,yvaluesub2,label=s)
    plt.axhline(y=0,color='black',linestyle='dashdot')
    plt.xlabel("Velocity(km/s)")
    plt.ylabel("Flux")
    plt.xlim(-20,5)
    plt.ylim(-0.1,0.1)
        
        
    
    plt.xlabel("Velocity (km/s)",fontsize=14)
    plt.ylabel("Flux (Jansky)",fontsize=14)
    plt.savefig("All",dpi=200,bbox_inches='tight')
    
    
    
    
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
        