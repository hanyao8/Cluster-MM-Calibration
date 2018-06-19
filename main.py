# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 15:49:41 2018

@author: Choon
"""

#main

#import io
import os
import numpy as np
#import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#import time
from datetime import datetime

import processing

###############USER INTERFACE#######################

#SETTING='general'
#SETTING='BD_angle'
#SETTING='mirrorratio'
#SETTING='themis_validation'
SETTING='peakness'

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

    
    params=[[data_start_time,data_end_time],[t_int,shift]]

    inst1 = processing.processing(df_arr,params,SETTING)
    

if SETTING=='themis_validation':
    import cdflib

    cdf_file1 = cdflib.CDF(os.getcwd()+"\\data\\thc_l2_fgm_20080702_v01.cdf") #02/07/2008 FGM level 1
    
    C_xy = 0.3
    C_B = 30*np.pi/180
    C_D = 30*np.pi/180
    C_phi = 20*np.pi/180
    
    t_int = 180
    shift = 10    
    
    data_start_time = matplotlib.dates.date2num(datetime.strptime('2008-07-02T14:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
    data_end_time = matplotlib.dates.date2num(datetime.strptime('2008-07-02T18:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))

    epoch1970 = 719163.0
    epoch20080702 = 733225.0
    
    fgs_Btot = cdf_file1.varget('thc_fgs_btotal')
    fgs_dsl = cdf_file1.varget('thc_fgs_dsl')
    fgs_t=cdf_file1.varget('thc_fgs_time')
    fgs_t2 = fgs_t/3600/24 + np.array([epoch1970]*len(fgs_t) ) 
    
    t = fgs_t2
    B_x = fgs_dsl[:,0]
    B_y = fgs_dsl[:,1]
    B_z = fgs_dsl[:,2]
    B_mag = fgs_Btot

    df_arr = [t,B_x,B_y,B_z,B_mag]
    params=[[data_start_time,data_end_time],[t_int,shift]]
    inst1 = processing.processing(df_arr,params,SETTING)
    
    Mmode_count = 0
    
    for i in range(0,len(inst1.Bxy_fluct_id)):
        #if Bxy_fluct_id[i]>C_xy:
        if abs(inst1.phi_PN16[i]) < C_phi:
            if abs(inst1.theta_B_PN16[i]) < C_B:
                if abs(inst1.theta_D_PN16[i]) < C_D:
                    Mmode_count+=1
    print(Mmode_count)

if SETTING=='BD_angle':
    #t_int_arr = [45,60,75,90,120,150]
    #t_int_arr = [15,30]
    #t_int_arr = [45,60,75]
    #t_int_arr = [90,120,150]
    t_int_arr=[300]
    #shift_arr=[15]
    #shift_arr = t_int_arr
    shift_arr = [10]*len(t_int_arr)

    f5=plt.figure()
    ax5 = f5.add_subplot(111)     

    data_start_time_str = '2006-03-01T10:30:00.000Z'
    data_end_time_str = '2006-03-01T11:29:00.000Z'
    
    data_start_time = matplotlib.dates.date2num(datetime.strptime(data_start_time_str,'%Y-%m-%dT%H:%M:%S.%fZ'))
    data_end_time = matplotlib.dates.date2num(datetime.strptime(data_end_time_str,'%Y-%m-%dT%H:%M:%S.%fZ'))

    var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']
    
    csv_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"

    csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_names)
    #csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_arr)
    csv_df.head()
    
    df_arr = csv_df.values
  
    for i in range(0,len(t_int_arr)):
        
        t_int = t_int_arr[i]
        shift = shift_arr[i]
    
        params=[[data_start_time,data_end_time],[t_int,shift]]
        inst1 = processing.processing(df_arr,params,SETTING)  
        
        sitv_midpoints = (inst1.subintervals[:,0]+inst1.subintervals[:,1])/2/3600/24
        
        ax5.plot_date(sitv_midpoints,inst1.Bx_angle,\
                       fmt='-',linewidth=1.0,label=r"$t_{si}$=%.1f,$t_{sh}$=%.1f"%(t_int,shift))
            
            
        ax5.plot_date(sitv_midpoints,\
                   np.array(len(sitv_midpoints)*[30]),fmt='-',linewidth=1.0)
        
        
    ax5.set_ylabel("B_D_angle (degs)")
    ax5.set_title("Time Series of Angle between MV dir. and B-field dir.")
    ax5.set_xlabel("Time")
    ax5.legend()
            
            
    plt.show()
            

if SETTING=='mirrorratio':
    from scipy.fftpack import fft
    
    t_int = 180
    shift = 60
    
    var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']
    
    csv_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"
    csv_df = pd.read_csv(os.getcwd()+"\\data\\"+  csv_file_name + ".csv",names=var_names)
    #csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_arr)
    csv_df.head()
    
    df_arr = csv_df.values

    params=[[data_start_time,data_end_time],[t_int,shift]]

    inst1 = processing.processing(df_arr,params,SETTING)
    
    
    
if SETTING=='peakness':

    
    t_int=300
    shift = 60
    
    data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T10:30:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
    data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T11:29:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
    
    var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']
    
    csv_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"

    csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_names)
    csv_df.head()
    
    df_arr = csv_df.values    
    #data of all variables within the interval defined above (data_start_time and data_end_time)

    
    params=[[data_start_time,data_end_time],[t_int,shift]]

    inst1 = processing.processing(df_arr,params,SETTING)

    inst1.b_wav
    
    
    
##########################PLOTTING ROUTINE###########################
    
if __name__=="__main__":
    if SETTING=='general':
        
        sitv_midpoints = (inst1.subintervals[:,0]+inst1.subintervals[:,1])/2/3600/24
        
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
        
        ax1.plot_date(inst1.t_days,inst1.B_mag,fmt='-',linewidth=1.0)
        ax1.set_title("B-field magnitude time series")
        ax1.set_xlabel("Time")
        ax1.set_ylabel(r"$B_{mag}$ (nT)")
        
        ax21.plot(inst1.t_days,inst1.B_x,linewidth=1.0)
        ax21.plot_date(sitv_midpoints,inst1.Bx_sitv_mean,fmt='-',linewidth=1.0)
        ax21.set_title("Cartesian B-field time series")
        ax21.set_ylabel(r"$B_{x}$ (nT)")
        ax22.plot(inst1.t_days,inst1.B_y,linewidth=1.0)
        ax22.plot_date(sitv_midpoints,inst1.By_sitv_mean,fmt='-',linewidth=1.0)
        ax22.set_ylabel(r"$B_{y}$ (nT)")
        ax23.plot(inst1.t_days,inst1.B_z,linewidth=1.0)
        ax23.plot_date(sitv_midpoints,inst1.Bz_sitv_mean,fmt='-',linewidth=1.0)
        ax23.set_ylabel(r"$B_{z}$ (nT)")
        ax23.set_xlabel("Time")
        
        
        ax31.plot_date(inst1.t_days,inst1.B_x/inst1.B_mag,fmt='-',linewidth=1.0)
        ax31.set_ylabel(r"$B_{x}$ (nT)")
        ax31.set_title("Normalised Cartesian B-field time series")
        ax32.plot_date(inst1.t_days,inst1.B_y/inst1.B_mag,fmt='-',linewidth=1.0)
        ax32.set_ylabel(r"$B_{y}$ (nT)")
        ax33.plot_date(inst1.t_days,inst1.B_z/inst1.B_mag,fmt='-',linewidth=1.0)
        ax33.set_ylabel(r"$B_{z}$ (nT)")
        ax33.set_xlabel("Time")
        
        
        
        ax40.plot_date(inst1.t_days,inst1.B_mag,fmt='-',linewidth=1.0)
        ax41.plot_date(inst1.t_days,np.array(inst1.theta_B)*180/np.pi,fmt='-',linewidth=1.0)
        ax42.plot_date(inst1.t_days,np.array(inst1.phi_B)*180/np.pi,fmt='-',linewidth=1.0)
        ax40.set_title("Sph. Polar B-field time series")
        ax40.set_ylabel(r"$B_{mag}$ (nT)")
        ax41.set_ylabel(r"$\theta$ (degs)")
        ax42.set_ylabel(r"$\phi$ (degs)")
        ax42.set_xlabel("Time")
        
        ax5.plot_date(sitv_midpoints,inst1.B_x1_angle,fmt='-',linewidth=1.0)
        ax5.plot_date(sitv_midpoints,np.array(len(sitv_midpoints)*[30]),fmt='-',linewidth=1.0)
        ax5.set_ylabel("B_D_angle (degs)")
        ax5.set_title("Angle between MV and B $t_{si}$=%d,$t_{sh}$=%d"%(t_int,shift))
        ax5.set_xlabel("Time")
        ax5.axvline(datetime(2006, 3, 1,11),alpha=0.5,color='orange')
        ax5.axvline(datetime(2006, 3, 1,11, 8),alpha=0.5,color='orange')
        
        ax61.plot_date(sitv_midpoints,np.array(inst1.theta_D)*180/np.pi,fmt='-',linewidth=1.0,label="MV")
        ax62.plot_date(sitv_midpoints,np.array(inst1.phi_D)*180/np.pi,fmt='-',linewidth=1.0)  
        ax61.set_title(r"Maximum Variance Direction $t_{si}$=%d,$t_{sh}$=%d"%(t_int,shift))
        ax61.set_ylabel(r"$\theta_{D}$ (degs)")
        ax62.set_ylabel(r"$\phi_{D}$ (degs)")
        ax61.plot_date(inst1.t_days,np.array(inst1.theta_B)*180/np.pi,fmt='-',linewidth=1.0,label="B")
        ax62.plot_date(inst1.t_days,np.array(inst1.phi_B)*180/np.pi,fmt='-',linewidth=1.0)    
        ax61.axvline(datetime(2006, 3, 1,11),alpha=0.6)
        ax61.axvline(datetime(2006, 3, 1,11, 8),alpha=0.6)
        ax62.axvline(datetime(2006, 3, 1,11),alpha=0.6)
        ax62.axvline(datetime(2006, 3, 1,11, 8),alpha=0.6)
        ax61.legend()
        
        plt.show()                




    if SETTING=='themis_validation':
        
        sitv_midpoints = (inst1.subintervals[:,0]+inst1.subintervals[:,1])/2/3600/24
        
        f1=plt.figure()
        f2=plt.figure()
        f3=plt.figure()
        f4=plt.figure()
        f5=plt.figure()
        f6=plt.figure()
        f7=plt.figure()
        
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
        ax71 = f7.add_subplot(211)
        ax72 = f7.add_subplot(212)    
        
        ax1.plot_date(inst1.t_days,B_mag,fmt='-',linewidth=1.0)
        ax1.set_title("B-field magnitude time series")
        ax1.set_xlabel("Time")
        ax1.set_ylabel(r"$B_{mag}$ (nT)")
        
        ax21.plot(inst1.t_days,B_x,linewidth=1.0)
        ax21.plot_date(sitv_midpoints,inst1.Bx_sitv_mean,fmt='-',linewidth=1.0)
        ax21.set_title("Cartesian B-field time series")
        ax21.set_ylabel(r"$B_{x}$ (nT)")
        ax22.plot(inst1.t_days,inst1.B_y,linewidth=1.0)
        ax22.plot_date(sitv_midpoints,inst1.By_sitv_mean,fmt='-',linewidth=1.0)
        ax22.set_ylabel(r"$B_{y}$ (nT)")
        ax23.plot(inst1.t_days,inst1.B_z,linewidth=1.0)
        ax23.plot_date(sitv_midpoints,inst1.Bz_sitv_mean,fmt='-',linewidth=1.0)
        ax23.set_ylabel(r"$B_{z}$ (nT)")
        ax23.set_xlabel("Time")
        
        
        ax31.plot_date(inst1.t_days,inst1.B_x/inst1.B_mag,fmt='-',linewidth=1.0)
        ax31.set_ylabel(r"$B_{x}$ (nT)")
        ax31.set_title("Normalised Cartesian B-field time series")
        ax32.plot_date(inst1.t_days,inst1.B_y/inst1.B_mag,fmt='-',linewidth=1.0)
        ax32.set_ylabel(r"$B_{y}$ (nT)")
        ax33.plot_date(inst1.t_days,inst1.B_z/inst1.B_mag,fmt='-',linewidth=1.0)
        ax33.set_ylabel(r"$B_{z}$ (nT)")
        ax33.set_xlabel("Time")
        
        
        """
        ax40.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)
        ax41.plot_date(t_days,np.array(theta)*180/np.pi,fmt='-',linewidth=1.0)
        ax42.plot_date(t_days,np.array(phi)*180/np.pi,fmt='-',linewidth=1.0)
        ax40.set_title("Sph. Polar B-field time series")
        ax40.set_ylabel(r"$B_{mag}$ (nT)")
        ax41.set_ylabel(r"$\theta$ (degs)")
        ax42.set_ylabel(r"$\phi$ (degs)")
        ax42.set_xlabel("Time")
        
        ax5.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,B_x1_angle,fmt='-',linewidth=1.0)
        ax5.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,np.array(len(subintervals[:,0])*[30]),fmt='-',linewidth=1.0)
        ax5.set_ylabel("B_D_angle (degs)")
        ax5.set_title("Time Series of Angle between MV dir. and B-field dir.")
        ax5.set_xlabel("Time")
        """
        
        ax61.plot_date(sitv_midpoints,np.abs(np.array(inst1.phi_PN16))*180/np.pi,fmt='-',linewidth=1.0)
        ax62.plot_date(sitv_midpoints,np.array(inst1.theta_D_PN16)*180/np.pi,fmt='-',linewidth=1.0)
        ax61.plot_date(sitv_midpoints,[20]*len(inst1.phi_PN16),fmt='-',linewidth=1.0)
        ax62.plot_date(sitv_midpoints,[30]*len(inst1.theta_D_PN16),fmt='-',linewidth=1.0)
        ax62.plot_date(sitv_midpoints,[-30]*len(inst1.theta_D_PN16),fmt='-',linewidth=1.0)
        ax61.set_ylabel(r"$\phi$")
        ax62.set_ylabel(r"$\theta_{D}$")
        
    
        ax71.plot_date(sitv_midpoints,np.array(inst1.theta_B_PN16)*180/np.pi,fmt='-',linewidth=1.0)
        ax71.plot_date(sitv_midpoints,[30]*len(inst1.theta_B_PN16),fmt='-',linewidth=1.0)
        ax71.plot_date(sitv_midpoints,[-30]*len(inst1.theta_B_PN16),fmt='-',linewidth=1.0)    
        ax72.plot_date(sitv_midpoints,inst1.Bxy_fluct_id,fmt='-',linewidth=1.0)    
        ax71.set_ylabel(r"$\theta_{B}$")
        ax72.set_ylabel(r"$\frac{\delta B_{xy}}{B_{xy}}$")
        
        plt.show()

    
    if SETTING=='mirrorratio':
        sitv_midpoints=((inst1.subintervals[:,0]+inst1.subintervals[:,1])/2/3600/24)
        
        f1=plt.figure()
        ax1=f1.add_subplot(111)
        f2=plt.figure()
        ax2=f2.add_subplot(111)
        
        a=np.array(len(sitv_midpoints)*[1.5])
        b=np.array(len(sitv_midpoints)*[0.3])
        
        
        mirrorindex=[]
        
        for i in range(len(sitv_midpoints)):
            rules=[inst1.lam1_lam2[i]>=1.5,inst1.lam3_lam2[i]>=0.3,inst1.B_x1_angle[i]<=30]
            if all(rules):
                mirrorindex.append(i)
                
        bf=fft(inst1.B_mag)
        ax1.plot(np.abs(bf[0:int(len(bf/2))]))
        
        ax2.plot_date(inst1.t_days,inst1.B_mag,fmt='-',linewidth=1.0)                
                
    if SETTING=='peakness':
        timeav=((inst1.subintervals[:,0]+inst1.subintervals[:,1])/2/3600/24)
        
        fig6=plt.figure()
        plt.plot_date(timeav,inst1.peakness,fmt='-',linewidth=1.0)
        plt.xlabel("Time")
        plt.ylabel("Peakness")
    



                                                                                     