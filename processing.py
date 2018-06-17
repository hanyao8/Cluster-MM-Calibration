# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 15:57:17 2018

@author: Choon
"""

#processing

import io
import os
import numpy as np
import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import time
from datetime import datetime


########################## FUNCTIONS #####################################

def varlist(data):
#Calculate the mean of the magnetic field components and return the variance matrix as a (3,3) array

    Bxm=data[0].mean()
    Bym=data[1].mean()
    Bzm=data[2].mean()
    Bxsqm=(data[0]**2).mean()
    Bysqm=(data[1]**2).mean()
    Bzsqm=(data[2]**2).mean()
    Bxym=(data[0]*data[1]).mean()
    Bxzm=(data[0]*data[2]).mean()
    Byzm=(data[1]*data[2]).mean()
    
    Varx= Bxsqm-Bxm**2
    Varxy=Bxym-Bxm*Bym
    Varxz=Bxzm-Bxm*Bzm
    Vary=Bysqm-Bym**2
    Varyz=Byzm-Bym*Bzm
    Varz=Bzsqm-Bzm**2
    var=[Varx,Varxy, Varxz,Varxy,Vary,Varyz,Varxz,Varyz,Varz]
    var=np.array(var)
    var=var.reshape((3,3))
    return var


def xyz2polar(B_x,B_y,B_z):
    B_mag = []
    theta = []
    phi = []
    for i in range(0,len(B_x)):
        B_mag.append( np.sqrt((B_x[i])**2 + (B_y[i])**2 + (B_z[i])**2 ) )
        theta.append( (np.arctan2( [np.sqrt((B_x[i])**2 + (B_y[i])**2 )],B_z[i] ))[0])
        phi.append( (np.arctan2( [B_y[i]],[B_x[i] ]) )[0] )
        #theta.append(np.arccos( np.dot(B_vec[i],[0,0,1]) ))
        #phi.append(np.arccos( np.dot( ([B_x[i],B_y[i],0]/B_xy[i]),[1,0,0]) ))
    return(np.array([B_mag,theta,phi]))

########################## FUNCTIONS #####################################



class processing:
    def __init__(self,df_arr,params):

        t = df_arr[:,0]
        B_x = df_arr[:,2]
        B_y = df_arr[:,3]
        B_z = df_arr[:,4]
        B_mag = df_arr[:,5]

        data_start_time = params[0][0]
        data_end_time = params[0][1]
        
        t_datetime = []
        for i in range(0,len(t)):
            strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
            t_datetime.append(strpdtime)
        t_days = matplotlib.dates.date2num(t_datetime)
        #print(t_days)
        t_secs= t_days*24*3600

        data_start_index = np.argmax(t_days>data_start_time)
        data_end_index = np.argmax(t_days>data_end_time)
        
        t = t[data_start_index:data_end_index]
        B_x = B_x[data_start_index:data_end_index]
        B_y = B_y[data_start_index:data_end_index]
        B_z = B_z[data_start_index:data_end_index]
        B_mag = B_mag[data_start_index:data_end_index]

        #t_datetime = []
        #for i in range(0,len(t)):
        #    strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
        #    t_datetime.append(strpdtime)
        #t_days = matplotlib.dates.date2num(t_datetime)
        #t_secs= t_days*24*3600
        
        t_days = t_days[data_start_index:data_end_index]
        t_secs = t_secs[data_start_index:data_end_index]

        
        B_xy = np.empty(len(t))
        for i in range(0,len(t)):
            B_xy[i] = np.sqrt((B_x[i])**2 + (B_y[i])**2 )

        B_vec = []
        for i in range(0,len(t)):
            B_vec.append(np.array([B_x[i],B_y[i],B_z[i]])/B_mag[i])    
    
        B_polar = xyz2polar(B_x,B_y,B_z)
        
        theta_B = B_polar[1]
        phi_B = B_polar[2]




        
        subintervals = np.empty(( int((t_secs[-1]-t_secs[0]-t_int+shift)/shift) ,2))
        for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
            subintervals[i][0] = t_secs[0] + i*shift
            subintervals[i][1] = t_secs[0] + i*shift + t_int

Bxy_sitv = []
Bx_sitv = []
By_sitv = []
Bz_sitv = []

Bx_sitv_mean = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))
By_sitv_mean = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))
Bz_sitv_mean = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))

Bxy_sitv_min = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))
Bxy_sitv_max = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))
Bxy_sitv_mean = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))

for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
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



Bx_angle = []
theta_D = []
phi_D = []
for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
    data = np.array([Bx_sitv[i],By_sitv[i],Bz_sitv[i]])
    M = varlist(data)
    eigen = np.linalg.eig(varlist(data))   
    
    lam1_index = np.argmax(eigen[0])
    
    x1 = eigen[1][:,lam1_index]
    x1_xy = np.sqrt(x1[0]**2+x1[1]**2)
    
    B_dir = np.array([Bx_sitv_mean[i],By_sitv_mean[i],Bz_sitv_mean[i]])
    B_dir /= np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2)
    B_dir_xy = np.sqrt(B_dir[0]**2+B_dir[1]**2)
    
    Bx_angle.append( np.arccos(np.dot(x1,B_dir)) * 180/np.pi )

    theta_D.append( (np.arctan2( [x1_xy],[x1[2]] ))[0])
    phi_D.append( (np.arctan2( [x1[1]],[x1[0] ]) )[0] )


