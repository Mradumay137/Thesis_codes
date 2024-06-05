import numpy as np
import re
import csv

timedata=np.loadtxt("//DATA/CETUS_7/sad043/mar829/P974/data/S37_data/time_data.txt", skiprows=6,dtype=str)
timecol=timedata[:,13]
source=timedata[:,6]
y=0
remove=[]
for x in range(0,len(timecol)):
    if source[x]=='S37':
        y=y+1
    if float(timecol[x])<20: 
        remove.append(y)        
rawfileS=np.loadtxt("//DATA/CETUS_7/sad043/mar829/P974/data/S37_data/S37.txt",dtype=str)
rawfileR=np.loadtxt("//DATA/CETUS_7/sad043/mar829/P974/data/S37_data/R2.txt",dtype=str)
rawS=rawfileS[:]
rawR=rawfileR[:]

with open("//DATA/CETUS_7/sad043/mar829/P974/data/S37_data/S37.txt") as f_input:
    text1 = [l.replace(".", "_") for l in f_input]
with open("//DATA/CETUS_7/sad043/mar829/P974/data/S37_data/R2.txt") as f_input:
    text2 = [l.replace(".", "_") for l in f_input]    
data1 = np.loadtxt(text1,delimiter="_",dtype=str)
data2 = np.loadtxt(text2,delimiter="_",dtype=str)

Sdate=data1[:,1]
Stime=data1[:,2]

Rdate=data2[:,1]
Rtime=data2[:,2]
optfile=[]

for i in range(0,len(Sdate)):
    k=1
    for j in range(0,len(Rdate)):
        datediff=int(Rdate[j])-int(Sdate[i])
        if datediff==0:
            diff=abs(int(Rtime[j])-int(Stime[i]))
            if diff<12000 and k<2 and i!=remove:
                k=k+1
                optfile.append(rawR[j])

rawS=np.asarray(rawS)
optfile=np.asarray(optfile)
if len(remove)>0:
    rawSdel=np.delete(rawS,remove)

with open('Matched_pairs.txt', 'w') as wf:
     writer=csv.writer(wf,delimiter=' ')
     if len(remove)>0:
         Column1=rawSdel
         Column2=optfile
     else:
         Column1=rawS
         Column2=optfile
     writer.writerows(zip(Column1,Column2))
with open('Matched_pairs.txt','r') as inf:  
    data =inf.read()
    corr=data.replace("^M"," ")
with open('Matched_pair_1.txt','w') as of:
    of.write(corr)
    
