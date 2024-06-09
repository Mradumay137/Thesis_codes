#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 21:28:49 2024

@author: maverick
"""
from scipy.optimize import curve_fit
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from PIL import Image, UnidentifiedImageError
plt.figure(figsize=(5,10))
for i  in range(1,3):
    data=np.loadtxt('Position.txt',skiprows=4)
    
    Rah=data[:,0]
    Ram=data[:,1]
    Ras=data[:,2]
    Rar=data[:,3]
    Decd=data[:,4]
    decam=data[:,5]
    decas=data[:,6]
    
    RA=(Rah + Ram/60 + Ras/3600)*15
    DEC=Decd-decam/60-decas/3600
    
    c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
    
    xcoord=[coord.l.deg for coord in c.galactic]
    ycoord=[coord.b.deg for coord in c.galactic]
    
    image_file = get_pkg_data_filename('highlat_iris100.fits')
    hdu = fits.open(image_file)[0]
    wcs = WCS(hdu.header)
    
    
    
    if i==1:
        ax1=plt.subplot(211,projection=wcs) 
    if i==2:
        ax2=plt.subplot(212,projection=wcs)     
    if i == 1:
        ax = ax1
    elif i == 2:
        ax = ax2
    
    marker_size_arcminutes = 12
    
    colordata=np.loadtxt('Output.txt')
    c1=colordata[:,0]
    c2=colordata[:,1]
    c3=colordata[:,2]
    c4=colordata[:,3]
    
    
    cdiff=[]
    cav=[]
    cperc=[]
    for gamma in range(0,len(c1)):
        cd=abs(abs(c1[gamma]-c2[gamma]))
        if c1[gamma]<0.85*10**14 or c2[gamma]<0.85*10**14:
            cperc.append(0)
        elif cd < abs(c3[gamma]):
            cperc.append(0)
        else:
            ca=(c1[gamma]+c2[gamma])/2
            cper=100*cd/ca
            cperc.append(cper)
            print(gamma)
    
    
    
    if i==1:
        ax.imshow(hdu.data, origin='lower') 
        ax.grid(color='white', ls='solid')
        plt.title("1667",loc='left')
        color=colordata[:,0]
        size=colordata[:,2]
        sizearr=[]
        percent1=[]
        file1=[]
        for beta in range(0,len(size)):
            siz=size[beta]/((10**11))
            """
            if color[beta]<1:
                sizearr.append(0)
            elif color[beta]>1:
            """    
            sizearr.append(siz)
            
            percentage1=(size[beta])*100/abs(color[beta])
            """
            if color[beta]<50 or size[beta]<50:
                percentage1=0  
            """
            percent1.append(percentage1)
            file1.append(beta)   
        s=ax.scatter(xcoord,ycoord,transform=ax.get_transform('world'), s=24,
               edgecolor='white', c=color,cmap='jet', vmin=0, vmax=2.8*10**14)
    if i==2:
        ax.imshow(hdu.data, origin='lower') 
        ax.grid(color='white', ls='solid')
        plt.title("1665",loc='left')
        color=colordata[:,1]
        size=colordata[:,3]
        
        sizearr=[]
        percent2=[]
        file2=[]
        for beta in range(0,len(size)):
            siz=size[beta]/((10**11))
            """
            if color[beta]<1:
                sizearr.append(0)
            elif color[beta]>1:
                """
            sizearr.append(siz)
            
            percentage2=(size[beta])*100/abs(color[beta])
            """
            if color[beta]<1 or size[beta]<1:
                percentage2=0    
            """
            if beta==23:
                percentage2=0
            percent2.append(percentage2)
            file2.append(beta)
             
        s=ax.scatter(xcoord,ycoord,transform=ax.get_transform('world'), s=24,
                   edgecolor='white', c=color,cmap='jet', vmin=0, vmax=2.8*10**14)
    def number_to_superscript(number):
        if number < 0:
            sign = "¯"
            number = abs(2)
        else:
            sign = " "
        superscript_digits = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
        return "{}{}".format(sign, str(2).translate(superscript_digits))
    
    # Example usage
    result = number_to_superscript(-2)
    plt.colorbar(s, ax=ax1 if i == 1 else ax2, label='OH column density (cm{})'.format(result))
plt.savefig("Poster",dpi=300)
plt.show()

#plt.figure(figsize=(15,15))

#plt.plot(file2,percent2,label='1665')
#xticks=np.arange(1,47,1) 
 
#plt.yticks(np.arange(0, 3000, step=100))  # Set label locations.

#plt.xticks(file2,xticks)
#plt.axhline(y=30)



#plt.plot(file1,percent1,label='1667') 
pixel_scale_deg = hdu.header['CDELT2']

# Convert degrees to arcseconds if needed
pixel_scale_arcsec = pixel_scale_deg * 3600

print(f'Pixel scale: {pixel_scale_arcsec:.2f} arcseconds per pixel')
#plt.legend()
"""
def func0(x,m,inter):
    return m*x+inter
plt.figure(figsize=(5,5))
plt.scatter(c1,c2)
plt.errorbar(c1,c2, yerr=c3,fmt="o")
plt.errorbar(c1,c2, xerr=c4,fmt="o")
x=np.linspace(0,3*10**14)
y=np.linspace(0,3*10**14)
plt.plot(x,y)
"""
"""
popt,pcov=curve_fit(func0,mean,size,p0=(2*10**14,2*10**13))
print(popt)

plt.figure(figsize=(5,5))
plt.scatter(mean, size)
plt.plot(mean,func0(mean,*popt),'k')
plt.xlabel('offset_from_mean')
plt.ylabel('error in integral')
"""
        