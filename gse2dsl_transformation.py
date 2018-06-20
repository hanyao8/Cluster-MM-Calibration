# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 15:02:33 2018

@author: Choon
"""

#gse2dsl transformation

import os
import numpy as np
#import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#import time
from datetime import datetime

SETTING='not_themis_validation'
#SETTING='themis_validation'

gsedatavar_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']
gsedata_file_name = "C3_CP_FGM_5VPS__20060301_000000_20060302_000000_V140305"
gsedata_df = pd.read_csv(os.getcwd()+"\\data\\"+ gsedata_file_name + ".csv",names=gsedatavar_names)
gsedata_df.head()
gsedata_arr = gsedata_df.values   

spinaxisvar_names = ['Time','Lat','Long']
spinaxis_file_name = "C3_CP_AUX_SPIN_AXIS__20060301_000000_20060302_000000_V130205"
spinaxis_df = pd.read_csv(os.getcwd()+"\\data\\"+ spinaxis_file_name + ".csv",names=spinaxisvar_names)
spinaxis_df.head()
spinaxis_arr = spinaxis_df.values   

df_arr = gsedata_arr.copy()
t = df_arr[:,0]
Bx_gse = df_arr[:,2]
By_gse = df_arr[:,3]
Bz_gse = df_arr[:,4]
B_mag = df_arr[:,5]


t_sa = spinaxis_arr[:,0]
lat_gse = spinaxis_arr[:,1]
long_gse = spinaxis_arr[:,2]


if SETTING=='themis_validation':
    t_days = t
else:
    t_datetime = []
    for i in range(0,len(t)):
        strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
        t_datetime.append(strpdtime)
    t_days = matplotlib.dates.date2num(t_datetime)
t_secs= t_days*24*3600


t_sa_datetime = []
for i in range(0,len(t_sa)):
    #print(t_sa[i])
    strpdtime=datetime.strptime(t_sa[i],'%Y-%m-%dT%H:%M:%S.%fZ')
    t_sa_datetime.append(strpdtime)
t_sa_days = matplotlib.dates.date2num(t_sa_datetime)
t_sa_secs= t_sa_days*24*3600


#converting spinaxis orientation time series into 5VPS resolution
lat_gse_newres = []
long_gse_newres = []
for i in range(0,len(t_days)):
    j=np.argmax(t_sa_days>t_days[i])
    lat_gse_newres.append((lat_gse[j-1]+lat_gse[j])/2)
    long_gse_newres.append((lat_gse[j-1]+lat_gse[j])/2)    
    
theta_gse = np.array([np.pi/2]*len(lat_gse_newres)) - np.array(lat_gse_newres)/180*np.pi
phi_gse = np.array(long_gse_newres)/180*np.pi

#Bx_dsl = []
#By_dsl = []
#Bz_dsl = []
for i in range(0,len(t_days)):
    z_dsl = np.array(([np.sin(theta_gse[i])*np.cos(phi_gse[i]),\
                  np.sin(theta_gse[i])*np.sin(phi_gse[i]),\
                  np.cos(theta_gse[i])]))
    x_cr_z = np.cross([1,0,0],z_dsl)
    y_dsl = ( x_cr_z / np.sqrt(np.dot(x_cr_z,x_cr_z)) )
    x_dsl = ( np.cross(y_dsl,z_dsl) )
    
    T = np.array([x_dsl,y_dsl,z_dsl])
    S = np.transpose(T) #matrix of 'eigenvectors'
    
    B_gse = [Bx_gse[i],By_gse[i],Bz_gse[i]]
    B_dsl = np.matmul(T,B_gse)
    
    
    
    






