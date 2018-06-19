# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 22:57:41 2018

@author: Choon
"""


#BD angle analysis

import io
import os
import numpy as np
import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import time
from datetime import datetime

#t_int_arr = [45,60,75,90,120,150]
#t_int_arr = [15,30]
#t_int_arr = [45,60,75]
#t_int_arr = [90,120,150]
t_int_arr=[300]
#shift_arr=[15]
#shift_arr = t_int_arr
shift_arr = [10]*len(t_int_arr)
#res = 0.2 #5 vectors per second

"""
f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
f4=plt.figure()

ax1 = f1.add_subplot(111)
ax21 = f2.add_subplot(311)
ax22 = f2.add_subplot(312)
ax23 = f2.add_subplot(313)
ax31 = f3.add_subplot(311)
ax32 = f3.add_subplot(312)
ax33 = f3.add_subplot(313)
ax40 = f4.add_subplot(311)
ax41 = f4.add_subplot(312)
ax42 = f4.add_subplot(313)
"""

f5=plt.figure()
ax5 = f5.add_subplot(111)                   

for i in range(0,len(t_int_arr)):
    
    t_int = t_int_arr[i]
    shift = shift_arr[i]
    
    data_start_time_str = '2006-03-01T10:30:00.000Z'
    data_end_time_str = '2006-03-01T11:29:00.000Z'
    
    data_start_time = matplotlib.dates.date2num(datetime.strptime(data_start_time_str,'%Y-%m-%dT%H:%M:%S.%fZ'))
    data_end_time = matplotlib.dates.date2num(datetime.strptime(data_end_time_str,'%Y-%m-%dT%H:%M:%S.%fZ'))
    
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
    
    var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']
    
    csv_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"
    #csv_data=np.genfromtxt(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",delimiter=',')
    
    
    
    
    ####################################################################
    
    
    #with open(os.getcwd()+"\\data\\"+ csv_file_name + ".csv") as csvfile:
    #    readCSV = csv.reader(csvfile,delimiter=',')
    
    csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_names)
    #csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_arr)
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
    
    
    data_start_index = np.argmax(t_days>data_start_time)
    data_end_index = np.argmax(t_days>data_end_time)
    t = t[data_start_index:data_end_index]
    B_x = B_x[data_start_index:data_end_index]
    B_y = B_y[data_start_index:data_end_index]
    B_z = B_z[data_start_index:data_end_index]
    B_mag = B_mag[data_start_index:data_end_index]
    
    t_datetime = []
    for i in range(0,len(t)):
        strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
        t_datetime.append(strpdtime)
    t_days = matplotlib.dates.date2num(t_datetime)
    
    t_secs= t_days*24*3600
    
    
    B_xy = np.empty(len(t))
    for i in range(0,len(t)):
        B_xy[i] = np.sqrt((B_x[i])**2 + (B_y[i])**2 )
    
    B_vec = []
    theta = []
    phi = []
    for i in range(0,len(t)):
        B_vec.append(np.array([B_x[i],B_y[i],B_z[i]])/B_mag[i])
        theta.append( (np.arctan2( [B_xy[i]],B_z[i] ))[0])
        phi.append( (np.arctan2( [B_y[i]],[B_x[i] ]) )[0] )
        #theta.append(np.arccos( np.dot(B_vec[i],[0,0,1]) ))
        #phi.append(np.arccos( np.dot( ([B_x[i],B_y[i],0]/B_xy[i]),[1,0,0]) ))
        
        
        
    
    
    
    ##################################################################################
    #t0_itv = t_secs[0]
    
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
    ##################################################################################
    """
    B_sitv=[Bx_sitv,By_sitv,Bz_sitv]
    
    M=np.empty((3,3))
    for i in range(0,3):
        for j in range(i,3):
            M[i][j]=M[j][i]=(np.cov(B_sitv[i],B_sitv[j]))[0][1]
    
    """
    ##################################################################################
    
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
    ##################################################################################
    
    Bx_angle = []
    for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
        data = np.array([Bx_sitv[i],By_sitv[i],Bz_sitv[i]])
        M = varlist(data)
        eigen = np.linalg.eig(varlist(data))   
        
        lam1_index = np.argmax(eigen[0])
        
        x1 = eigen[1][:,lam1_index]
        
        B_dir = np.array([Bx_sitv_mean[i],By_sitv_mean[i],Bz_sitv_mean[i]])
        B_dir /= np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2)
        
        Bx_angle.append( np.arccos(np.dot(x1,B_dir)) * 180/np.pi )
    
    ##################################################################################
    
    
    
    
        
        
     
    ##################################################################################

if __name__=="__main__":

    """
    ax1.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)
    ax1.set_title("B-field magnitude time series")
    ax1.set_xlabel("Time")
    ax1.set_ylabel(r"$B_{mag}$ (nT)")
    
    ax21.plot(t_days,B_x,linewidth=1.0)
    ax21.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,Bx_sitv_mean,fmt='-',linewidth=1.0)
    ax21.set_title("Cartesian B-field time series")
    ax21.set_ylabel(r"$B_{x}$ (nT)")
    ax22.plot(t_days,B_y,linewidth=1.0)
    ax22.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,By_sitv_mean,fmt='-',linewidth=1.0)
    ax22.set_ylabel(r"$B_{y}$ (nT)")
    ax23.plot(t_days,B_z,linewidth=1.0)
    ax23.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,Bz_sitv_mean,fmt='-',linewidth=1.0)
    ax23.set_ylabel(r"$B_{z}$ (nT)")
    ax23.set_xlabel("Time")
    
    
    ax31.plot_date(t_days,B_x/B_mag,fmt='-',linewidth=1.0)
    ax31.set_ylabel(r"$B_{x}$ (nT)")
    ax31.set_title("Normalised Cartesian B-field time series")
    ax32.plot_date(t_days,B_y/B_mag,fmt='-',linewidth=1.0)
    ax32.set_ylabel(r"$B_{y}$ (nT)")
    ax33.plot_date(t_days,B_z/B_mag,fmt='-',linewidth=1.0)
    ax33.set_ylabel(r"$B_{z}$ (nT)")
    ax33.set_xlabel("Time")
    
    
    
    ax40.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)
    ax41.plot_date(t_days,np.array(theta)*180/np.pi,fmt='-',linewidth=1.0)
    ax42.plot_date(t_days,np.array(phi)*180/np.pi,fmt='-',linewidth=1.0)
    ax40.set_title("Sph. Polar B-field time series")
    ax40.set_ylabel(r"$B_{mag}$ (nT)")
    ax41.set_ylabel(r"$\theta$ (degs)")
    ax42.set_ylabel(r"$\phi$ (degs)")
    ax42.set_xlabel("Time")
    """
    
    ax5.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,Bx_angle,\
                   fmt='-',linewidth=1.0,label=r"$t_{si}$=%.1f,$t_{sh}$=%.1f"%(t_int,shift))
        
        
    ax5.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,\
               np.array(len(subintervals[:,0])*[30]),fmt='-',linewidth=1.0)
    
    
ax5.set_ylabel("B_D_angle (degs)")
ax5.set_title("Time Series of Angle between MV dir. and B-field dir.")
ax5.set_xlabel("Time")
ax5.legend()
        
        
plt.show()
        
    

