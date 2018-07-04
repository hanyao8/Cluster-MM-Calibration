#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 16:24:34 2018

@author: jia_qu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 10:27:04 2018

@author: jia_qu
"""

import os
import numpy as np
#import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#import time
from datetime import datetime
from scipy import stats
import gaussian
from scipy.integrate import quad

from sklearn.neighbors.kde import KernelDensity
import cdflib
ARTIFICIAL_OFFSET = False
PLOT= True
sitv_in_sheath_frac_analysis = True

artificial_Bz_offset = 5.0 #nT (in DSL coordinates)

#csv_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"
cdf_file1 = cdflib.CDF("/Users/jia_qu/Desktop/thc_l2_fgm_20080702_v01.cdf") #02/07/2008 FGM level 1

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
"""
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
"""



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
#B_mag = df_arr[:,5]

#!!!!!!!!!!!!! artificial offset set here:
if ARTIFICIAL_OFFSET:
    B_z = np.array(len(B_z)*[artificial_Bz_offset]) + np.array(B_z)

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
t_days = t2

t_secs= t_days*24*3600

"""
self.t_days=t_days
self.t_secs=t_secs
"""
#print(t_days[0])
#print(self.t_days[0])


B_xy = np.empty(len(B_x))
for i in range(0,len(B_x)):
    B_xy[i] = np.sqrt((B_x[i])**2 + (B_y[i])**2 )

B_vec = []
for i in range(0,len(B_x)):
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


B_x1_angle = []
theta_D = []
phi_D = []
lam1_lam2 = []
lam3_lam2 = []
lam2_lam1=[]

Bxy_fluct_PN16 = (np.array(Bxy_sitv_max)-np.array(Bxy_sitv_min))/Bxy_sitv_mean
phi_PN16 = []
theta_B_PN16 = []
theta_D_PN16 = []

O_z = []
dO=[]

d_thetaD=[]
d_g=10**(-4)
d_bn=10**(-11)
d_b=[]
bxymirror=[]
bzmirror=[]
d_thetaB=[]
theta_mb=[]
theta_md=[]
test=[]
B_magn=[]
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
    lam2_lam1.append(lam2/lam1)
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
    theta_D.append( (np.arctan2([x1_xy],[x1[2]]))[0])
    phi_D.append( (np.arctan2( [x1[1]],[x1[0] ]) )[0] )
    B_magn.append(np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2))
    
    #if SETTING=='mmode_search':
    B_x1_angle.append( np.arccos(abs(np.dot(x1,B_dir)))  )
    phi_PN16.append( np.arccos(abs(np.dot([x1[0],x1[1]],[B_dir[0],B_dir[1]]) \
                                     /np.sqrt(x1[0]**2+x1[1]**2)/np.sqrt(B_dir[0]**2+B_dir[1]**2) ))) 
    
    theta_B_PN16.append( np.arctan2(B_dir[2],B_dir_xy) )
    theta_D_PN16.append( np.arctan2(x1[2],x1_xy) )
    rules=[Bxy_fluct_PN16[i]>C_xy,abs(phi_PN16[i]) < C_phi,abs(theta_D_PN16[i]) < C_D]
    if Bxy_fluct_PN16[i]>C_xy:
        if abs(phi_PN16[i]) < C_phi:
            if abs(theta_B_PN16[i]) < C_B:
                if abs(theta_D_PN16[i]) < C_D:
                     d_thetaD.append(np.arctan(lam2_lam1[i]))
                     d_b.append(abs(B_magn[i])*d_g+d_bn)
                     #d_b.append(abs(Bz_sitv_mean[i])*d_g+d_bn)
                     bxymirror.append(Bxy_sitv_mean[i])
                     bzmirror.append(Bz_sitv_mean[i])
                     theta_mb.append(theta_B_PN16[i])
                     theta_md.append(theta_D_PN16[i])
                     test.append(Bxy_sitv_mean[i]*(np.tan(theta_B_PN16[i])-np.tan(theta_D_PN16[i])))
                     O_z.append(Bz_sitv_mean[i]-x1[2]/x1_xy*Bxy_sitv_mean[i])
                     
                     """
                     if Bz_sitv_mean[i] > 0:
                        O_z.append(+abs(Bz_sitv_mean[i]) - abs(x1[2])/x1_xy*Bxy_sitv_mean[i] )
                     elif Bz_sitv_mean[i] < 0:
                        O_z.append(-abs(Bz_sitv_mean[i]) + abs(x1[2])/x1_xy*Bxy_sitv_mean[i] )
 
                    """
error=[]
for i in range(len(d_b)):
    d_thetaB.append((d_b[i]/(1+(bzmirror[i]/bxymirror[i])**2))*np.sqrt((1/bxymirror[i])**2+(bzmirror[i]/bxymirror[i]**2)**2))   
    error.append(np.sqrt(((np.tan(theta_mb[i])-np.tan(theta_md[i]))*d_b[i])**2+((bxymirror[i]*d_thetaB[i])/(np.cos(theta_mb[i]))**2)**2+((bxymirror[i]*d_thetaD[i])/(np.cos(theta_md[i]))**2)**2))                 
test=np.array(test)

O_z = np.array(O_z)
O_z_T = np.transpose(np.array([O_z]))

O_z_T = np.transpose(np.array([O_z]))
kde= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(O_z_T)
log_pdf = kde.score_samples(O_z_T)
plt.figure()
plt.scatter(O_z,log_pdf,s=1.0,color='orange',label="0nT artificial offset")
plt.figure()
x=np.linspace(0,30,100)
y2=x
y=-1*x
plt.plot(error,O_z,'+')
plt.plot(x,y)
plt.plot(x,y2)