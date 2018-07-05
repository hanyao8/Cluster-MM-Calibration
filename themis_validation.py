# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 10:02:57 2018

@author: Choon
"""

#mva validation themis

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
import cdflib


cdf_file1 = cdflib.CDF(os.getcwd()+"\\data\\thc_l2_fgm_20080702_v01.cdf") #02/07/2008 FGM level 1

#print(cdf_file1.cdf_info())

C_xy = 0.3
C_B = 30*np.pi/180
C_D = 30*np.pi/180
C_phi = 20*np.pi/180

t_int = 180
shift = 10
#res = 0.2 #5 vectors per second



data_start_time = matplotlib.dates.date2num(datetime.strptime('2008-07-02T14:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
data_end_time = matplotlib.dates.date2num(datetime.strptime('2008-07-02T18:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))


####################################################################


"""
csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_names)
csv_df.head()

df_arr = csv_df.values
"""
"""
t = df_arr[:,0]
B_x = df_arr[:,2]
B_y = df_arr[:,3]
B_z = df_arr[:,4]
B_mag = df_arr[:,5]
"""

epoch1970 = 719163.0
epoch20080702 = 733225.0

fgs_Btot = cdf_file1.varget('thc_fgs_btotal')
fgs_dsl = cdf_file1.varget('thc_fgs_dsl')
fgs_t=cdf_file1.varget('thc_fgs_time')
#fgs_t2 = (fgs_t - np.array([fgs_t[0]]*len(fgs_t))) /3600/24
fgs_t2 = fgs_t/3600/24 + np.array([epoch1970]*len(fgs_t) ) 


#matplotlib.dates.num2date(1214957297.0989337/3600/24 + epoch1970)

#t = fgs_t2 + np.array([epoch20080702]*len(fgs_t))
t = fgs_t2
B_x = fgs_dsl[:,0]
B_y = fgs_dsl[:,1]
B_z = fgs_dsl[:,2]
B_mag = fgs_Btot


"""
t_datetime = []
for i in range(0,len(t)):
    strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
    t_datetime.append(strpdtime)
"""
#t_days = matplotlib.dates.date2num(t_datetime)
t_days = t

#t_secs= t_days*24*3600


data_start_index = np.argmax(t_days>data_start_time)
data_end_index = np.argmax(t_days>data_end_time)
t2 = t[data_start_index:data_end_index]
B_x = B_x[data_start_index:data_end_index]
B_y = B_y[data_start_index:data_end_index]
B_z = B_z[data_start_index:data_end_index]
B_mag = B_mag[data_start_index:data_end_index]

"""
t_datetime = []
for i in range(0,len(t)):
    strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
    t_datetime.append(strpdtime)
"""
#t_days = matplotlib.dates.date2num(t_datetime)
t_days = t2

t_secs= t_days*24*3600


B_xy = np.empty(len(t_days))
for i in range(0,len(t_days)):
    B_xy[i] = np.sqrt((B_x[i])**2 + (B_y[i])**2 )

B_vec = []
theta = []
phi = []
for i in range(0,len(t_days)):
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
    #Bxy_sitv_mean[i] = np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2)
    
    
##################################################################################

#B_sitv=[Bx_sitv,By_sitv,Bz_sitv]

#M=np.empty((3,3))
#for i in range(0,3):
#    for j in range(i,3):
#        M[i][j]=M[j][i]=(np.cov(B_sitv[i],B_sitv[j]))[0][1]


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

Bxy_fluct_id = (Bxy_sitv_max-Bxy_sitv_min)/Bxy_sitv_mean
phi_id = []
theta_B_id = []
theta_D_id = []
#phi_id = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))
#= [np.arccos(abs(np.dot([Bx_sitv_mean[i],By_sitv_mean[i]],[Bx_sitv_mean[i],By_sitv_mean[i]]))) for i in range(0,len(Bx_sitv_mean))]
#theta_B_id = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))
#theta_D_id = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))

Bx_angle = []
for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
    data = np.array([Bx_sitv[i],By_sitv[i],Bz_sitv[i]])
    M = varlist(data)
    eigen = np.linalg.eig(varlist(data))   
    
    lam1_index = np.argmax(eigen[0])
    #lam1_index = np.argmin(eigen[0])
    
    x1 = eigen[1][:,lam1_index]
    x1_xy = np.sqrt(x1[0]**2+x1[1]**2)
    
    
    
    if i==0:
        print(x1)
        x1 = +1*np.array([-0.27334004 , 0.89995506 , 0.33965588])
        #x1 = +x1
        #x1 = -x1
        x1_prev = x1.copy()
        
    else:
        if np.dot(x1_prev,x1) < 0:
            x1 = -x1
    
    x1_prev = x1.copy()
    
    B_dir = np.array([Bx_sitv_mean[i],By_sitv_mean[i],Bz_sitv_mean[i]])
    B_dir /= np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2)
    B_dir_xy = np.sqrt(B_dir[0]**2+B_dir[1]**2)
    
    Bx_angle.append( np.arccos(abs(np.dot(x1,B_dir))) * 180/np.pi )

    phi_id.append( np.arccos(np.dot([x1[0],x1[1]],[B_dir[0],B_dir[1]]) \
                                     /np.sqrt(x1[0]**2+x1[1]**2)/np.sqrt(B_dir[0]**2+B_dir[1]**2) ) )
    
    """
    phi_id_val = ( np.arctan2(x1[1],x1[0]) - np.arctan2(B_dir[1],B_dir[0]) ) 
    if phi_id_val > np.pi:
        phi_id.append(phi_id_val-2*np.pi)
    elif phi_id_val < -np.pi:
        phi_id.append(phi_id_val+2*np.pi)
    else:
        phi_id.append(phi_id_val)
    """
    theta_B_id.append( np.arctan2(B_dir[2],B_dir_xy) )
    theta_D_id.append( np.arctan2(x1[2],x1_xy) )

##################################################################################





Mmode_count = 0

for i in range(0,len(Bxy_fluct_id)):
    #if Bxy_fluct_id[i]>C_xy:
    if abs(phi_id[i]) < C_phi:
        if abs(theta_B_id[i]) < C_B:
            if abs(theta_D_id[i]) < C_D:
                Mmode_count+=1


    
    
 
##################################################################################
if __name__=="__main__":
    
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
    
    
    """
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
    ax5.set_title("Time Series of Angle between MV dir. and B-field dir.")
    ax5.set_xlabel("Time")
    """
    
    ax61.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,np.abs(np.array(phi_id))*180/np.pi,fmt='-',linewidth=1.0)
    ax62.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,np.array(theta_D_id)*180/np.pi,fmt='-',linewidth=1.0)
    ax61.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,[20]*len(phi_id),fmt='-',linewidth=1.0)
    ax62.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,[30]*len(theta_D_id),fmt='-',linewidth=1.0)
    ax62.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,[-30]*len(theta_D_id),fmt='-',linewidth=1.0)
    ax61.set_ylabel(r"$\phi$")
    ax62.set_ylabel(r"$\theta_{D}$")
    

    ax71.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,np.array(theta_B_id)*180/np.pi,fmt='-',linewidth=1.0)
    ax71.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,[30]*len(theta_B_id),fmt='-',linewidth=1.0)
    ax71.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,[-30]*len(theta_B_id),fmt='-',linewidth=1.0)    
    ax72.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,Bxy_fluct_id,fmt='-',linewidth=1.0)    
    ax71.set_ylabel(r"$\theta_{B}$")
    ax72.set_ylabel(r"$\frac{\delta B_{xy}}{B_{xy}}$")
    
    plt.show()
    


"""
"""