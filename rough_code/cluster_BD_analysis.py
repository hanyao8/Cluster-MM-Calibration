# -*- coding: utf-8 -*-
"""
Created on Sun Jun 10 16:36:47 2018

@author: Choon
"""

#cluster_BD_analysis

import io
import os
import numpy as np
import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import time
from datetime import datetime

t_int = 180
overlap = 10
res = 0.2 #5 vectors per second

var_arr = ['time_tags__C1_CP_FGM_5VPS',
 'half_interval__C1_CP_FGM_5VPS',
 'B_vec_x_gse__C1_CP_FGM_5VPS',
 'B_vec_y_gse__C1_CP_FGM_5VPS',
 'B_vec_z_gse__C1_CP_FGM_5VPS',
 'B_mag__C1_CP_FGM_5VPS',
 'sc_pos_x_gse__C1_CP_FGM_5VPS',
 'sc_pos_y_gse__C1_CP_FGM_5VPS',
 'sc_pos_z_gse__C1_CP_FGM_5VPS', 
 'range__C1_CP_FGM_5VPS',
 'tm__C1_CP_FGM_5VPS']

csv_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"
#csv_data=np.genfromtxt(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",delimiter=',')

#with open(os.getcwd()+"\\data\\"+ csv_file_name + ".csv") as csvfile:
#    readCSV = csv.reader(csvfile,delimiter=',')

csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_arr)
csv_df.head()

df_arr = csv_df.values

t = df_arr[:,0]
B_x = df_arr[:,2]
B_y = df_arr[:,3]
B_z = df_arr[:,4]
B_mag = df_arr[:,5]


t_datetime = []
for i in range(0,len(t)):
    strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
    t_datetime.append(strpdtime)
t_days = matplotlib.dates.date2num(t_datetime)

t_secs= t_days*24*3600

B_xy = np.empty(len(t_secs))
for i in range(0,len(t_secs)):
    B_xy[i] = np.sqrt((B_x[i])**2 + (B_y[i])**2 )





##################################################################################
t0_itv = t_secs[0]

subintervals = np.empty(( int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap) ,2))
for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap)):
    subintervals[i][0] = t0_itv + i*overlap
    subintervals[i][1] = t0_itv + i*overlap + t_int

Bxy_sitv = []
Bx_sitv = []
By_sitv = []
Bz_sitv = []

Bx_sitv_mean = np.empty(int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap))
By_sitv_mean = np.empty(int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap))
Bz_sitv_mean = np.empty(int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap))

Bxy_sitv_min = np.empty(int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap))
Bxy_sitv_max = np.empty(int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap))
Bxy_sitv_mean = np.empty(int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap))

for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap)):
    Bx_sitv.append(B_x[np.argmax(t_secs>subintervals[i][0]):np.argmax(t_secs>subintervals[i][1])])
    By_sitv.append(B_y[np.argmax(t_secs>subintervals[i][0]):np.argmax(t_secs>subintervals[i][1])])
    Bz_sitv.append(B_z[np.argmax(t_secs>subintervals[i][0]):np.argmax(t_secs>subintervals[i][1])])
    
    Bxy_sitv.append(B_xy[np.argmax(t_secs>subintervals[i][0]):np.argmax(t_secs>subintervals[i][1])])
    
    
    Bx_sitv_mean[i] = np.mean(Bx_sitv[i])
    By_sitv_mean[i] = np.mean(By_sitv[i])
    Bz_sitv_mean[i] = np.mean(Bz_sitv[i])    
    
    
    Bxy_sitv_min[i] = min(Bxy_sitv[i])
    Bxy_sitv_max[i] = max(Bxy_sitv[i])
Bxy_sitv_mean[i] = np.mean(Bxy_sitv[i])
##################################################################################

B_sitv=[Bx_sitv,By_sitv,Bz_sitv]

M=np.empty((3,3))
for i in range(0,3):
    for j in range(i,3):
        M[i][j]=M[j][i]=(np.cov(B_sitv[i],B_sitv[j]))[0][1]


        


##################################################################################


f1=plt.figure()
f2=plt.figure()
ax1 = f1.add_subplot(111)
ax2 = f2.add_subplot(111)

ax1.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)

ax2.plot(t_secs,B_z,linewidth=1.0)
ax2.plot((subintervals[:,0]+subintervals[:,1])/2,Bz_sitv_mean,linewidth=1.0)

plt.show()




