# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 15:49:41 2018

@author: Choon
"""

#main

import io
import os
import numpy as np
import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import time
from datetime import datetime

import processing

###############USER INTERFACE#######################

SETTING='general'
SETTING='BD_angle'


if SETTING=='general':
    t_int = 180
    #t_int = 300
    #shift=t_int
    shift = 180
    #res = 0.2 #5 vectors per second
    
    data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T10:30:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
    data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T11:29:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
    
    var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']
    
    csv_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"

    csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_names)
    csv_df.head()
    
    df_arr = csv_df.values    
    #data of all variables within the interval defined above (data_start_time and data_end_time)

    
    params=[[data_start_time,data_end_time]]

    inst1 = processing.processing(df_arr,params)
    





if __name__=="__main__":
    if SETTING=='general':
        
        f1=plt.figure()
        f2=plt.figure()
        f3=plt.figure()
        f4=plt.figure()
        f5=plt.figure()
        f6=plt.figure()
        
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
        ax5 = f5.add_subplot(111)
        ax61 = f6.add_subplot(211)
        ax62 = f6.add_subplot(212)
        
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
        
        ax5.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,Bx_angle,fmt='-',linewidth=1.0)
        ax5.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,np.array(len(subintervals[:,0])*[30]),fmt='-',linewidth=1.0)
        ax5.set_ylabel("B_D_angle (degs)")
        ax5.set_title("Angle between MV and B $t_{si}$=%d,$t_{sh}$=%d"%(t_int,shift))
        ax5.set_xlabel("Time")
        ax5.axvline(datetime(2006, 3, 1,11),alpha=0.5,color='orange')
        ax5.axvline(datetime(2006, 3, 1,11, 8),alpha=0.5,color='orange')
        
        ax61.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,np.array(theta_D)*180/np.pi,fmt='-',linewidth=1.0,label="MV")
        ax62.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,np.array(phi_D)*180/np.pi,fmt='-',linewidth=1.0)  
        ax61.set_title(r"Maximum Variance Direction $t_{si}$=%d,$t_{sh}$=%d"%(t_int,shift))
        ax61.set_ylabel(r"$\theta_{D}$ (degs)")
        ax62.set_ylabel(r"$\phi_{D}$ (degs)")
        ax61.plot_date(t_days,np.array(theta)*180/np.pi,fmt='-',linewidth=1.0,label="B")
        ax62.plot_date(t_days,np.array(phi)*180/np.pi,fmt='-',linewidth=1.0)    
        ax61.axvline(datetime(2006, 3, 1,11),alpha=0.6)
        ax61.axvline(datetime(2006, 3, 1,11, 8),alpha=0.6)
        ax62.axvline(datetime(2006, 3, 1,11),alpha=0.6)
        ax62.axvline(datetime(2006, 3, 1,11, 8),alpha=0.6)
        ax61.legend()
        
        plt.show()                                                                                                     