#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:09:03 2024

@author: maverick
"""

import numpy as np 
import matplotlib.pyplot as plt

data=np.loadtxt('int_intensity.txt')

col1=data[:,0]

col2=np.array([1,2,3,4])

col3=np.array([6,7,8,10])

col4=data[:,1]
# Chunk size
chunk_size = 4

# List to hold the chunks
chunked_lists1 = []
chunked_lists2 = []
# Using a loop to iterate over the original list and create chunks
for i in range(0, len(col1), chunk_size):
    chunk = col1[i:i + chunk_size]
    chunked_lists1.append(chunk)
for i in range(0, len(col4), chunk_size):
    chunk = col4[i:i + chunk_size]
    chunked_lists2.append(chunk)    

I = [chunk[:4] for chunk in chunked_lists1]
err= [chunk[:4] for chunk in chunked_lists2]

bar_width = 1

colors = ['blue', 'green', 'orange', 'red']
# Plotting each sublist
for i, sublist in enumerate(I):
    # Create a new figure for each subplot
    plt.scatter(col2,sublist,color=colors[i])
    plt.errorbar(col2, sublist, yerr=err[i], label=f'S{col3[i]}',color=colors[i],barsabove=True,capsize=10)  # Plot the sublist
plt.xlabel('Method')  # Set x-axis label
plt.ylabel('Integrated intensity')  
plt.xticks(range(1,len(I)+1))
plt.title('Jo and George plot')  # Set title
plt.legend()  # Add legend
plt.grid(True)  # Add grid
plt.show()  # Show the plot
