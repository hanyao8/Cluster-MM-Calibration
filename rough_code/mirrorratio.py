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
from scipy.fftpack import fft
t_int = 180
overlap = 60
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

var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']

csv_file_name = "200603"
#csv_data=np.genfromtxt(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",delimiter=',')

#with open(os.getcwd()+"\\data\\"+ csv_file_name + ".csv") as csvfile:
#    readCSV = csv.reader(csvfile,delimiter=',')

csv_df = pd.read_csv(os.getcwd()+"//" +  csv_file_name + ".csv",names=var_names)
#csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_arr)
csv_df.head()

df_arr = csv_df.values

t = df_arr[1:,0]
B_x = df_arr[1:,2]
B_y = df_arr[1:,3]
B_z = df_arr[1:,4]
B_mag = df_arr[1:,5]


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
ratiomaxmin=[]
ratiomaxmid=[]
ratiomidmin=[]
lmax=[]
lmin=[]
lmid=[]

for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+overlap)/overlap)):
    data = np.array([Bx_sitv[i],By_sitv[i],Bz_sitv[i]])
    M = varlist(data)
    eigen = np.linalg.eig(varlist(data))   
    lam1_index = np.argmax(eigen[0])
    #recheck the code for degenerate situations
    eigmax=np.max(eigen[0])
    lmax.append(eigmax)
    eigmin=np.min(eigen[0])
    lmin.append(eigmin)
    for i in range(len(eigen[0])):
        if eigen[0][i]!=eigmax or eigmin:
            eigmid=eigen[0][i]
    lmid.append(eigmid)
    ratiomaxmin.append(eigmax/eigmin)
    ratiomaxmid.append(eigmax/eigmid)
    ratiomidmin.append(eigmin/eigmid)
    x1 = eigen[1][:,lam1_index]
    
    B_dir = np.array([Bx_sitv_mean[i],By_sitv_mean[i],Bz_sitv_mean[i]])
    B_dir /= np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2)
    
    Bx_angle.append( np.arccos((np.dot(-x1,B_dir))) * 180/np.pi )
    
##################################################################################

    
    
    
 
##################################################################################


"""
f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
ax1 = f1.add_subplot(111)
ax2 = f2.add_subplot(111)
ax3 = f3.add_subplot(111)

ax1.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)

ax2.plot(t_secs,B_z,linewidth=1.0)
ax2.plot((subintervals[:,0]+subintervals[:,1])/2,Bz_sitv_mean,linewidth=1.0)

ax3.plot((subintervals[:,0]+subintervals[:,1])/2,Bx_angle,linewidth=1.0)

plt.show()
"""
a=np.array(len(ratiomaxmid)*[1.5])
b=np.array(len(ratiomaxmid)*[0.3])
##################################################################################
"""
plt.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,ratiomaxmid,fmt='-',label="max-mid eigenratio")
plt.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,ratiomaxmin,fmt='-',label="max-min eigenratio")
plt.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,ratiomidmin,fmt='-',label="min-mid eigenratio")
plt.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,a,fmt='-')
plt.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,b,fmt='-')

plt.title("Eigenvalues ratio vs time 2006 March")
plt.xlabel("interval/s")
plt.legend()
fig1=plt.figure()
plt.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)
f2=plt.figure()
plt.plot(t_secs,B_z,linewidth=1.0)
plt.plot((subintervals[:,0]+subintervals[:,1])/2,Bz_sitv_mean,linewidth=1.0)
fig3=plt.figure()
plt.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,lmax,fmt='-',label="eigmax")
plt.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,lmin,fmt='-',label="eigmin")
plt.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,lmid,fmt='-',label="eigmid")
plt.title("Eigenvalues of MVA vs time 2006 March")
plt.ylabel("variance/$nT^2$")
plt.xlabel("interval/s")
plt.legend()
"""

timeav=((subintervals[:,0]+subintervals[:,1])/2/3600/24)
mirrorindex=[]

for i in range(len(ratiomaxmid)):
    rules=[ratiomaxmid[i]>=1.5,ratiomidmin[i]>=0.3,Bx_angle[i]<=30]
    if all(rules):
        mirrorindex.append(i)
"""
fig4=plt.figure()
plt.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)
for xc in mirrorindex:
    print(timeav[xc])
    plt.axvline(x=timeav[xc],color='k')

fig5=plt.figure()
plt.plot_date(timeav,Bx_angle,fmt='-')
c=np.array(len(ratiomaxmid)*[30])
plt.plot_date((subintervals[:,0]+subintervals[:,1])/2/3600/24,c,fmt='-')
"""


bf=fft(B_mag)
plt.plot(np.abs(bf[0:int(len(bf/2))]))

plt.figure()
plt.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)