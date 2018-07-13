#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 11:30:54 2018

@author: jia_qu
"""

#seath identification
import os

import numpy as np
#import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#import time
from datetime import datetime
import io
import pywt
from scipy.integrate import quad


t_int = 300
shift = 300

data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-02-01T00:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-02-28T23:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))


var_names=['time_tags__C3_CP_CIS-HIA_ONBOARD_MOMENTS',
 'delta_time__C3_CP_CIS-HIA_ONBOARD_MOMENTS',
 'sensitivity__C3_CP_CIS-HIA_ONBOARD_MOMENTS',
 'cis_mode__C3_CP_CIS-HIA_ONBOARD_MOMENTS',
 'density__C3_CP_CIS-HIA_ONBOARD_MOMENTS',
 'velocity_isr2x',
 'velocity_isr2y'
 'velocity_isr2z'
 'velocity_gsex',
 'velocity_gsey',
 'velocity_gsez',
 'temperature__C3_CP_CIS-HIA_ONBOARD_MOMENTS',
 'temp_par__C3_CP_CIS-HIA_ONBOARD_MOMENTS',
 'temp_perp__C3_CP_CIS-HIA_ONBOARD_MOMENTS',
 'pressure__C3_CP_CIS-HIA_ONBOARD_MOMENTS',
 'pressure_tensor__C3_CP_CIS-HIA_ONBOARD_MOMENTS']

csv_file_name = "C3_CP_CIS-HIA_ONBOARD_MOMENTS__20060201_000000_20060301_000000_V161018"
fgm_file="C3_CP_FGM_SPIN__20060201_000000_20060301_000000_V140305"
csv_df = pd.read_csv(os.getcwd()+"//" +  csv_file_name + ".csv")
fgm_df=pd.read_csv(os.getcwd()+"//" +  fgm_file+ ".csv")

df_arr = csv_df.values    
fgm_arr=fgm_df.values
iond=df_arr[:,4]
t = df_arr[:,0]
t_datetime = []
for i in range(0,len(t)):
    strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
    t_datetime.append(strpdtime)
t_days = matplotlib.dates.date2num(t_datetime)
t_secs= t_days*24*3600

data_start_index = np.argmax(t_days>data_start_time)
data_end_index = np.argmax(t_days>data_end_time)
     

t = t[data_start_index:data_end_index]
t_days = t_days[data_start_index:data_end_index]
t_secs = t_secs[data_start_index:data_end_index]

subintervals = np.empty(( int((t_secs[-1]-t_secs[0]-t_int+shift)/shift) ,2))
for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
    subintervals[i][0] = t_secs[0] + i*shift
    subintervals[i][1] = t_secs[0] + i*shift + t_int
  
gapadj_subintervals = []
for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
    si_start = np.argmax(t_secs>subintervals[i][0])
    si_end = np.argmax(t_secs>subintervals[i][1])

iondensitysitv=[]
for i in range(0,len(subintervals)):
    si_start = np.argmax(t_secs>subintervals[i][0])
    si_end = np.argmax(t_secs>subintervals[i][1])  
    
    iondensitysitv.append(np.mean(iond[si_start:si_end]))

ionindex=[]
if iondensitysitv[0]>2*7.1 and iondensitysitv[1]>2*7.1:
    ionindex.append(0)
elif iondensitysitv[-1]>2*7.1 and iondensitysitv[-2]>2*7.1:
    ionindex.append(len(iondensitysitv)-1)
for i in range(1,len(iondensitysitv)-1):
    rules=[iondensitysitv[i-1]>2*7.1,iondensitysitv[i]>2*7.1,iondensitysitv[i+1]>2*7.1]
    if all(rules):
        ionindex.append(i)

sitv_midpoints_days = (subintervals[:,0]+subintervals[:,1])/2/3600/24
sheathtime=[]
for i in range(0,len(ionindex)):
    sheathtime.append(matplotlib.dates.num2date(sitv_midpoints_days[ionindex[i]]))
    #Mmode_sitv_times_PN16_days.append(sitv_midpoints_days[Mmode_indices_PN16[i]])
        

f1=plt.figure()
ax1 = f1.add_subplot(111)
a=np.array([2*7.1]*len(iond))
ax1.plot_date(t_datetime,iond,fmt='-')
ax1.plot_date(t_datetime,a,fmt='-',color='r')
#ax1.set_title("B-field magnitude time series")
ax1.set_xlabel("Time")
ax1.set_ylabel("Ion density/$cm^{-3}$")


for i in range(0,len(ionindex)):          
    ax1.axvline(sheathtime[i],alpha=0.2,color='orange',linewidth=1.0)
    
subtime=[]
for i in ionindex:
    
    subtime.append([subintervals[i][0]/(24*60*60),subintervals[i][1]/(24*60*60)])
    
i=0
j=0 
newgfmfile=io.open(os.getcwd()+"//"+"newfgmfile"+".csv","a")

print(len(fgm_arr[:,0]))
while i <len(fgm_arr[:,0]):
    
    if i%10000==0:
        print(i)
        
    anything=matplotlib.dates.date2num(datetime.strptime(fgm_arr[:,0][i],'%Y-%m-%dT%H:%M:%S.%fZ'))
    if anything>subtime[j][0] and anything<subtime[j][1]:
        data_row_str=fgm_arr[i][0]+","+str(fgm_arr[i][1])+str(fgm_arr[i][2])+","+str(fgm_arr[i][3])+","+str(fgm_arr[i][4])
        newgfmfile.write(data_row_str+'\n')
        i+=1
    elif anything < subtime[j][0]:
        i+=1
    elif anything > subtime[j][1]:
        if j<len(ionindex)-1:
            j+=1
        elif j==len(ionindex)-1:
            break

 

