# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 16:30:14 2018

@author: Choon
"""

#calibration

#import io
import os
import numpy as np
#import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#import time
from datetime import datetime

import pywt
from scipy.integrate import quad

from sklearn.neighbors.kde import KernelDensity

#offset_analysis
#contains code from 'mmode_search' part of main.py
#dsl based data must be used for this to give valid results

################################################################################
    
######################## UI AND SETTINGS ##############################################

################################################################################

ARTIFICIAL_OFFSET = False
PLOT= True
sitv_in_sheath_frac_analysis = True

artificial_Bz_offset = 5.0 #nT (in DSL coordinates)

#csv_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"
csv_file_name = "DSL_C3_CP_FGM_5VPS__20060301_000000_20060302_000000_V140305"

data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T00:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T23:59:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))

t_int = 180
shift = 10

C_xy = 0.3
C_B = 30*np.pi/180
C_D = 30*np.pi/180
C_phi = 20*np.pi/180

C_lam12 = 1.5
C_lam32 = 0.3
C_MV_B = 30*np.pi/180

var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']

csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_names)
csv_df.head()

df_arr = csv_df.values


################################################################################
    
######################## FUNCTIONS ##############################################

################################################################################

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


################################################################################
    
######################## PROCESSING ##############################################

################################################################################


t = df_arr[:,0]
B_x = df_arr[:,2]
B_y = df_arr[:,3]
B_z = df_arr[:,4]
#B_mag = df_arr[:,5]

#!!!!!!!!!!!!! artificial offset set here:
if ARTIFICIAL_OFFSET:
    B_z = np.array(len(B_z)*[artificial_Bz_offset]) + np.array(B_z)

"""
if SETTING=='themis_validation':
    t_days = t

else:
"""    
    
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
#B_mag = B_mag[data_start_index:data_end_index]
B_x.astype(float)
B_y.astype(float)
B_z.astype(float)
B_mag = (B_x**2+B_y**2+B_z**2)**(0.5)

#print(t[0])
"""
self.t=t

self.B_x=B_x
self.B_y=B_y
self.B_z=B_z
self.B_mag=B_mag
"""

#t_datetime = []
#for i in range(0,len(t)):
#    strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
#    t_datetime.append(strpdtime)
#t_days = matplotlib.dates.date2num(t_datetime)
#t_secs= t_days*24*3600

t_days = t_days[data_start_index:data_end_index]
t_secs = t_secs[data_start_index:data_end_index]
"""
self.t_days=t_days
self.t_secs=t_secs
"""
#print(t_days[0])
#print(self.t_days[0])


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


Bx_sitv = []
By_sitv = []
Bz_sitv = []
Bmag_sitv=[]
Bxy_sitv = []

Bx_sitv_mean = []
By_sitv_mean = []
Bz_sitv_mean = []

Bxy_sitv_min = []
Bxy_sitv_max = []
Bxy_sitv_mean = []

gapadj_subintervals = []
for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
    si_start = np.argmax(t_secs>subintervals[i][0])
    si_end = np.argmax(t_secs>subintervals[i][1])

    if len(B_x[si_start:si_end])!=0:
        gapadj_subintervals.append(subintervals[i])
        
        Bx_sitv.append(B_x[si_start:si_end])
        By_sitv.append(B_y[si_start:si_end])
        Bz_sitv.append(B_z[si_start:si_end])
        Bmag_sitv.append(B_mag[si_start:si_end])
        Bxy_sitv.append(B_xy[si_start:si_end])
        
        Bx_sitv_mean.append(np.mean(B_x[si_start:si_end]))
        By_sitv_mean.append(np.mean(B_y[si_start:si_end]))
        Bz_sitv_mean.append(np.mean(B_z[si_start:si_end]))
        """
        if i%1000==0:
            print(si_start)
            print(si_end)
            print(np.mean(self.B_x[si_start:si_end]))
            #print(self.Bx_sitv_mean[i])
        """
        Bxy_sitv_min.append(min(Bxy_sitv[-1]))
        Bxy_sitv_max.append(max(Bxy_sitv[-1]))
        Bxy_sitv_mean.append(np.mean(Bxy_sitv[-1]))

Bx_sitv_mean = np.array(Bx_sitv_mean)
By_sitv_mean = np.array(By_sitv_mean)
Bz_sitv_mean = np.array(Bz_sitv_mean)             
        
subintervals = np.array(gapadj_subintervals.copy())
#self.subintervals = subintervals

"""
print(self.Bx_sitv_mean)
"""

B_x1_angle = []
theta_D = []
phi_D = []
lam1_lam2 = []
lam3_lam2 = []

Bxy_fluct_PN16 = (np.array(Bxy_sitv_max)-np.array(Bxy_sitv_min))/Bxy_sitv_mean
phi_PN16 = []
theta_B_PN16 = []
theta_D_PN16 = []

Mmode_count_PN16 = 0
Mmode_count_SLD08 = 0

Mmode_indices_PN16 = []
Mmode_indices_SLD08 = []
Mmode_indices_overlapping = []



O_z = []
#using PN16 criteria 

for i in range(0,len(subintervals)):
    data = np.array([Bx_sitv[i],By_sitv[i],Bz_sitv[i]])
    eigen = np.linalg.eig(varlist(data))   

    lam1_index = np.argmax(eigen[0])
    lam3_index = np.argmin(eigen[0])
    for eigenindex in [0,1,2]:
        if eigenindex!=lam1_index and eigenindex!=lam3_index:
            lam2_index=eigenindex
            break
        
    lam1 = eigen[0][lam1_index]
    lam2 = eigen[0][lam2_index]
    lam3 = eigen[0][lam3_index]
    
    #print(lam1,lam2,lam3)
    
    lam1_lam2.append(lam1/lam2)
    lam3_lam2.append(lam3/lam2)

    x1 = eigen[1][:,lam1_index]
    x1_xy = np.sqrt(x1[0]**2+x1[1]**2)

    x1_prev = [0.0,0.0,0.0] #just an initialisation to make python happy
    if i==0:
        #x1 = +1*np.array([-0.27334004 , 0.89995506 , 0.33965588])
        x1 = +x1
        #x1 = -x1
        #x1_prev = x1.copy()
        
    else:
        if np.dot(x1_prev,x1) < 0:
            x1 = -x1
    
    x1_prev = x1.copy()
    



    B_dir = np.array([Bx_sitv_mean[i],By_sitv_mean[i],Bz_sitv_mean[i]])
    B_dir /= np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2)
    B_dir_xy = np.sqrt(B_dir[0]**2+B_dir[1]**2)
    
    """
    if i%1000==0:            
        print(i)
        print(self.Bx_sitv_mean[i])
        print(self.Bx_sitv_mean)
        print(B_dir_xy)
        print(np.arctan2(B_dir[2],B_dir_xy))     
    """

    theta_D.append( (np.arctan2( [x1_xy],[x1[2]] ))[0])
    phi_D.append( (np.arctan2( [x1[1]],[x1[0] ]) )[0] )
    
    #if SETTING=='mmode_search':
    B_x1_angle.append( np.arccos(abs(np.dot(x1,B_dir)))  )
    phi_PN16.append( np.arccos(abs(np.dot([x1[0],x1[1]],[B_dir[0],B_dir[1]]) \
                                     /np.sqrt(x1[0]**2+x1[1]**2)/np.sqrt(B_dir[0]**2+B_dir[1]**2) ))) 
    """           
    else:
        B_x1_angle.append( np.arccos(np.dot(x1,B_dir))  )
        phi_PN16.append( np.arccos(np.dot([x1[0],x1[1]],[B_dir[0],B_dir[1]]) \
                                         /np.sqrt(x1[0]**2+x1[1]**2)/np.sqrt(B_dir[0]**2+B_dir[1]**2) ) )
    """
    #if i%1000==0:
        #print(B_dir[2])
        #print(B_dir_xy)
        #print(np.arctan2(B_dir[2],B_dir_xy))
    theta_B_PN16.append( np.arctan2(B_dir[2],B_dir_xy) )
    theta_D_PN16.append( np.arctan2(x1[2],x1_xy) )
    
    
    if Bxy_fluct_PN16[i]>C_xy:
        if abs(phi_PN16[i]) < C_phi:
            if abs(theta_B_PN16[i]) < C_B:
                if abs(theta_D_PN16[i]) < C_D:
                    Mmode_indices_PN16.append(i)
                    Mmode_count_PN16+=1
                    
                    if Bz_sitv_mean[i] > 0:
                        O_z.append(+abs(Bz_sitv_mean[i]) - abs(x1[2])/x1_xy*Bxy_sitv_mean[i] )
                    elif Bz_sitv_mean[i] < 0:
                        O_z.append(-abs(Bz_sitv_mean[i]) + abs(x1[2])/x1_xy*Bxy_sitv_mean[i] )
                    
    if lam1_lam2[i] > C_lam12:
        if lam3_lam2[i] > C_lam32:
            if B_x1_angle[i] < C_MV_B:
                Mmode_indices_SLD08.append(i)
                Mmode_count_SLD08+=1
                
                if len(Mmode_indices_PN16)>0 and len(Mmode_indices_SLD08)>0:
                    if Mmode_indices_PN16[-1]==Mmode_indices_SLD08[-1]:
                        Mmode_indices_overlapping.append(i)    

################################################################################
    
######################## PLOTTING ##############################################

################################################################################

O_z = np.array(O_z)




sitv_midpoints_days = (subintervals[:,0]+subintervals[:,1])/2/3600/24
Mmode_sitv_times_PN16 = []
Mmode_sitv_times_SLD08 = []
Mmode_sitv_times_PN16_days = []
Mmode_sitv_times_SLD08_days = []

for i in range(0,len(Mmode_indices_PN16)):
    Mmode_sitv_times_PN16.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_PN16[i]]))
    Mmode_sitv_times_PN16_days.append(sitv_midpoints_days[Mmode_indices_PN16[i]])
for i in range(0,len(Mmode_indices_SLD08)):
    Mmode_sitv_times_SLD08.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_SLD08[i]]))   
    Mmode_sitv_times_SLD08_days.append(sitv_midpoints_days[Mmode_indices_SLD08[i]])

Mmode_sitv_times_PN16_only = Mmode_sitv_times_PN16
Mmode_sitv_times_SLD08_only = Mmode_sitv_times_SLD08
Mmode_sitv_times_overlapping = [] 
Mmode_sitv_times_overlapping_days = [] 
for i in range(0,len(Mmode_indices_overlapping)):        
    Mmode_sitv_times_overlapping.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_overlapping[i]])) 
    Mmode_sitv_times_overlapping_days.append(sitv_midpoints_days[Mmode_indices_overlapping[i]])
    Mmode_sitv_times_PN16_only.remove(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_overlapping[i]]))
    Mmode_sitv_times_SLD08_only.remove(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_overlapping[i]]))
    
    
if sitv_in_sheath_frac_analysis:
    sheath_start_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T06:58:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
    sheath_end_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T11:40:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
    
    PN16_frac_in_sheath = (np.argmax(np.array(Mmode_sitv_times_PN16_days)>sheath_end_time) - np.argmax(np.array(Mmode_sitv_times_PN16_days)>sheath_start_time) )/Mmode_count_PN16
    SLD08_frac_in_sheath = (np.argmax(np.array(Mmode_sitv_times_SLD08_days)>sheath_end_time) - np.argmax(np.array(Mmode_sitv_times_SLD08_days)>sheath_start_time) )/Mmode_count_SLD08  
    overlapping_frac_in_sheath = (np.argmax(np.array(Mmode_sitv_times_overlapping_days)>sheath_end_time) - np.argmax(np.array(Mmode_sitv_times_overlapping_days)>sheath_start_time) )/len(Mmode_indices_overlapping)
    
    print("PN16:",Mmode_count_PN16)
    print("SLD08:",Mmode_count_SLD08)
    print("overlapping:",len(Mmode_indices_overlapping))
    print("PN16 only:",len(Mmode_sitv_times_PN16_only))
    print("SLD08 only:",len(Mmode_sitv_times_SLD08_only))
    
    print(PN16_frac_in_sheath)
    print(SLD08_frac_in_sheath)
    print(overlapping_frac_in_sheath)



#Obtaining the pdf using kde
O_z_T = np.transpose(np.array([O_z]))
kde= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(O_z_T)
log_pdf = kde.score_samples(O_z_T)

        










################################################################################
    
######################## PLOTTING ##############################################

################################################################################

if PLOT:
    
    
    sitv_midpoints = (subintervals[:,0]+subintervals[:,1])/2/3600/24
    
    f1=plt.figure()
    ax1 = f1.add_subplot(111)
    
    f2=plt.figure()
    ax21 = f2.add_subplot(311)
    ax22 = f2.add_subplot(312)
    ax23 = f2.add_subplot(313)
    
    f3=plt.figure()
    ax31 = f3.add_subplot(311)
    ax32 = f3.add_subplot(312)
    ax33 = f3.add_subplot(313)
    
    #f4=plt.figure()
    #ax40 = f4.add_subplot(311)
    #ax41 = f4.add_subplot(312)
    #ax42 = f4.add_subplot(313)
    
    f5=plt.figure()
    ax5 = f5.add_subplot(111)
    
    f6=plt.figure()
    ax61 = f6.add_subplot(211)
    ax62 = f6.add_subplot(212)
    
    f7=plt.figure()
    ax71 = f7.add_subplot(211)
    ax72 = f7.add_subplot(212)        
    
    f8=plt.figure()
    ax81 = f8.add_subplot(211)
    ax82 = f8.add_subplot(212)      

    f9=plt.figure()
    ax9 = f9.add_subplot(111)
    
    
    ax1.plot_date(t_days,B_mag,fmt='-',linewidth=1.0,color='black')
    ax1.set_title("B-field magnitude time series")
    ax1.set_xlabel("Time")
    ax1.set_ylabel(r"$B_{mag}$ (nT)")
    """
    for i in range(0,len(Mmode_sitv_times_PN16)):
        #i%10==0:            
        ax1.axvline(Mmode_sitv_times_PN16[i],alpha=0.4,color='orange',linewidth=1.0)
    for i in range(0,len(Mmode_sitv_times_SLD08)):
        ax1.axvline(Mmode_sitv_times_SLD08[i],alpha=0.3,color='green',linewidth=1.0)    
    ax1.axvline(Mmode_sitv_times_PN16[0],alpha=0.4,color='orange',linewidth=1.0,label='PN16')
    ax1.axvline(Mmode_sitv_times_SLD08[0],alpha=0.3,color='green',linewidth=1.0,label='SLD08')
    """
    for i in range(0,len(Mmode_sitv_times_PN16_only)):          
        ax1.axvline(Mmode_sitv_times_PN16_only[i],alpha=0.4,color='orange',linewidth=1.0)
    for i in range(0,len(Mmode_sitv_times_SLD08_only)):
        ax1.axvline(Mmode_sitv_times_SLD08_only[i],alpha=0.4,color='green',linewidth=1.0)    
    for i in range(0,len(Mmode_sitv_times_overlapping)):
        ax1.axvline(Mmode_sitv_times_overlapping[i],alpha=0.4,color='blue',linewidth=1.0)    
    ax1.axvline(Mmode_sitv_times_PN16[0],alpha=0.4,color='orange',linewidth=1.0,label='PN16 only')
    ax1.axvline(Mmode_sitv_times_SLD08[0],alpha=0.4,color='green',linewidth=1.0,label='SLD08(no[4]) only')    
    ax1.axvline(Mmode_sitv_times_overlapping[0],alpha=0.4,color='blue',linewidth=1.0,label='PN16&SLD08(no[4])')    
    ax1.legend()
    
    ax21.plot(t_days,B_x,linewidth=1.0)
    ax21.plot_date(sitv_midpoints,Bx_sitv_mean,fmt='-',linewidth=1.0)
    ax21.set_title("Cartesian B-field time series")
    ax21.set_ylabel(r"$B_{x}$ (nT)")
    ax22.plot(t_days,B_y,linewidth=1.0)
    ax22.plot_date(sitv_midpoints,By_sitv_mean,fmt='-',linewidth=1.0)
    ax22.set_ylabel(r"$B_{y}$ (nT)")
    ax23.plot(t_days,B_z,linewidth=1.0)
    ax23.plot_date(sitv_midpoints,Bz_sitv_mean,fmt='-',linewidth=1.0)
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
    
    """
    ax40.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)
    ax41.plot_date(t_days,np.array(theta)*180/np.pi,fmt='-',linewidth=1.0)
    ax42.plot_date(t_days,np.array(phi)*180/np.pi,fmt='-',linewidth=1.0)
    ax40.set_title("Sph. Polar B-field time series")
    ax40.set_ylabel(r"$B_{mag}$ (nT)")
    ax41.set_ylabel(r"$\theta$ (degs)")
    ax42.set_ylabel(r"$\phi$ (degs)")
    ax42.set_xlabel("Time")
    """
    ax5.plot_date(sitv_midpoints,np.array(B_x1_angle)*180/np.pi,fmt='-',linewidth=1.0)
    ax5.plot_date(sitv_midpoints,np.array(len(sitv_midpoints)*[30]),fmt='-',linewidth=1.0)
    ax5.set_ylabel("B_D_angle (degs)")
    ax5.set_title("Time Series of Angle between MV dir. and B-field dir.")
    ax5.set_xlabel("Time")
    
    
    ax61.plot_date(sitv_midpoints,np.abs(np.array(phi_PN16))*180/np.pi,fmt='-',linewidth=1.0)
    ax62.plot_date(sitv_midpoints,np.array(theta_D_PN16)*180/np.pi,fmt='-',linewidth=1.0)
    ax61.plot_date(sitv_midpoints,[20]*len(phi_PN16),fmt='-',linewidth=1.0)
    ax62.plot_date(sitv_midpoints,[30]*len(theta_D_PN16),fmt='-',linewidth=1.0)
    ax62.plot_date(sitv_midpoints,[-30]*len(theta_D_PN16),fmt='-',linewidth=1.0)
    ax61.set_ylabel(r"$\phi$")
    ax62.set_ylabel(r"$\theta_{D}$")
    
    
    ax71.plot_date(sitv_midpoints,np.array(theta_B_PN16)*180/np.pi,fmt='-',linewidth=1.0)
    ax71.plot_date(sitv_midpoints,[30]*len(theta_B_PN16),fmt='-',linewidth=1.0)
    ax71.plot_date(sitv_midpoints,[-30]*len(theta_B_PN16),fmt='-',linewidth=1.0)    
    ax72.plot_date(sitv_midpoints,Bxy_fluct_PN16,fmt='-',linewidth=1.0)    
    ax72.plot_date(sitv_midpoints,[C_xy]*len(Bxy_fluct_PN16),fmt='-',linewidth=1.0)          
    ax71.set_ylabel(r"$\theta_{B}$")
    ax72.set_ylabel(r"$\frac{\delta B_{xy}}{B_{xy}}$")
    
    ax81.plot_date(sitv_midpoints,np.array(lam1_lam2),fmt='-',linewidth=1.0)
    ax81.plot_date(sitv_midpoints,[C_lam12]*len(lam1_lam2),fmt='-',linewidth=1.0) 
    ax82.plot_date(sitv_midpoints,np.array(lam3_lam2),fmt='-',linewidth=1.0)    
    ax82.plot_date(sitv_midpoints,[C_lam32]*len(lam3_lam2),fmt='-',linewidth=1.0)   
    ax81.set_ylabel(r"$\frac{\lambda_{1}}{\lambda_{2}}$")
    ax82.set_ylabel(r"$\frac{\lambda_{3}}{\lambda_{2}}$")     
    
    #ax9.hist(O_z,bins=25,normed=True)
    if ARTIFICIAL_OFFSET:
        ax9.scatter(O_z,np.exp(log_pdf),s=1.0,color='orange',label="%.2fnT artificial offset"%(artificial_Bz_offset))
    else:
        ax9.scatter(O_z,np.exp(log_pdf),s=1.0,color='orange',label="0nT artificial offset")
    ax9.set_title("Prob. Density from KDE: C3 01/03/2006, %d PN16 intervals"%(Mmode_count_PN16))
    ax9.set_ylabel("Prob. Density")
    ax9.set_xlabel(r"$O_{z}$ (nT)")
    ax9.legend()
    
    plt.show()        






