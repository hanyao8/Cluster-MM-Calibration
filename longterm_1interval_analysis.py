# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 16:30:14 2018

@author: Choon
"""

#calibration

#import io
import os
import numpy as np
import math
#import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#import time
from datetime import datetime
import bandwidth

#import pywt
#from scipy.integrate import quad

from sklearn.neighbors.kde import KernelDensity
import mmcal_functions

#offset_analysis
#contains code from 'mmode_search' part of main.py
#dsl based data must be used for this to give valid results

################################################################################
    
######################## UI AND SETTINGS ##############################################

################################################################################

ARTIFICIAL_OFFSET = False
PLOT= True
sitv_in_sheath_frac_analysis = False
simulate_spinres = False

artificial_Bz_offset = 5.0 #nT (in DSL coordinates- hence in spin axis direction)
if not(ARTIFICIAL_OFFSET):
    artificial_Bz_offset = 0.0

#csv_file_name = "DSL_C3_CP_FGM_5VPS__20061103_000000_20061104_000000_V140305"
#csv_file_name = "DSL_C3_CP_FGM_5VPS__20060301_000000_20060302_000000_V140305"
#csv_file_name = "DSL_C3_CP_FGM_5VPS__20060302_000000_20060307_000000_V140305"
#csv_file_name = "DSL_C3_CP_FGM_5VPS__20060307_000000_20060312_000000_V140305"
csv_file_name = "DSL_THEMIS_C_FGM_SPINRES_2008Jul"

#data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T00:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T22:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T21:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T22:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T00:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T23:59:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))

yyyy_s = '2008'
mm_s = '07'
dd_s = '02'
hh24_s = '14'

yyyy_e = '2008'
mm_e = '07'
dd_e = '02'
hh24_e = '18'

plot_time_text = yyyy_s+'/'+mm_s+'/'+dd_s+' '+hh24_s+'~'+hh24_e

data_start_time = matplotlib.dates.date2num(datetime.strptime(yyyy_s+'-'+mm_s+'-'+dd_s+'T'+hh24_s+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
data_end_time = matplotlib.dates.date2num(datetime.strptime(yyyy_e+'-'+mm_e+'-'+dd_e+'T'+hh24_e+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))

sheath_start_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T06:58:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
sheath_end_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T11:40:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))

plot_time_text = "Feb 2006"

t_int = 180
shift = 10

C_xy = 0.3
C_B = 30*np.pi/180
C_D = 30*np.pi/180
C_phi = 20*np.pi/180

C_lam12 = 1.5
C_lam32 = 0.3
C_MV_B = 30*np.pi/180

#var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']
var_names = ['Time','Half Interval','Bx','By','Bz']

#csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_names)
csv_df = pd.read_csv(os.getcwd()+"\\"+ csv_file_name + ".csv",names=var_names)
csv_df.head()
df_arr = csv_df.values
print("Progress: Main DSL data file loaded")

################################################################################

######################## END OF UI #############################################

################################################################################

################################################################################

######################## FUNCTIONS #############################################

################################################################################

if csv_file_name[:3]!='DSL':
    print("Ensure B-field data is in DSL coordinates")



def fivevps_2_spinres(t_overall,B_x_overall,B_y_overall,B_z_overall):
    t_overall_new = []
    B_x_overall_new = []
    B_y_overall_new = []
    B_z_overall_new = []   
    for i in range(0,len(t_overall)):
        if i%20==0:
            t_overall_new.append(t_overall[i])
            B_x_overall_new.append(B_x_overall[i])
            B_y_overall_new.append(B_y_overall[i])    
            B_z_overall_new.append(B_z_overall[i])
    
    return( [t_overall_new,B_x_overall_new,B_y_overall_new,B_z_overall_new] )

################################################################################
    
######################## PROCESSING ##############################################

################################################################################

t = df_arr[:,0]
B_x = df_arr[:,2]
B_y = df_arr[:,3]
B_z = df_arr[:,4]
#B_mag = df_arr[:,5]

if simulate_spinres:
    data_overall_spinres = fivevps_2_spinres(t,B_x,B_y,B_z)
    t_overall = data_overall_spinres[0]
    B_x_overall = np.array(data_overall_spinres[1])
    B_y_overall = np.array(data_overall_spinres[2])
    B_z_overall = np.array(data_overall_spinres[3])        

#!!!!!!!!!!!!! artificial offset set here:

if ARTIFICIAL_OFFSET:
    B_z = np.array(len(B_z)*[artificial_Bz_offset]) + np.array(B_z)


O_z_overall_PN16 = []
O_z_overall_overlapping = []
error_O_z_PN16 = []
error_O_z_overlapping = []


modified_data = mmcal_functions.data_in_interval([t,B_x,B_y,B_z],[data_start_time,data_end_time])
t = modified_data[0][0]
t_days = modified_data[0][1]
t_secs = modified_data[0][2]
B_x = modified_data[1][0]
B_y = modified_data[1][1]
B_z = modified_data[1][2]
B_mag = modified_data[2]
B_xy =np.array( [(np.sqrt((B_x[i])**2 + (B_y[i])**2 ) ) for i in range(0,len(t))] )

B_vec =np.array( [(np.array([B_x[i],B_y[i],B_z[i]])/B_mag[i]) for i in range(0,len(t))]) 

B_polar = mmcal_functions.xyz2polar(B_x,B_y,B_z)

theta_B = B_polar[1]
phi_B = B_polar[2]


subintervals = np.empty(( int((t_secs[-1]-t_secs[0]-t_int+shift)/shift) ,2))
for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
    subintervals[i][0] = t_secs[0] + i*shift
    subintervals[i][1] = t_secs[0] + i*shift + t_int


sitvfied = mmcal_functions.sitvfy([t_secs,[B_x,B_y,B_z],B_mag],[t_int,shift])
sitv_midpoints_secs = sitvfied[0]
Bx_sitv = sitvfied[1][0]
By_sitv = sitvfied[1][1]
Bz_sitv = sitvfied[1][2]
Bx_sitv_mean = [np.mean(Bx_sitv[i]) for i in range(0,len(sitv_midpoints_secs))]
By_sitv_mean = [np.mean(By_sitv[i]) for i in range(0,len(sitv_midpoints_secs))]
Bz_sitv_mean = [np.mean(Bz_sitv[i]) for i in range(0,len(sitv_midpoints_secs))]
Bmag_sitv_mean = sitvfied[2]
Bxy_sitv_min = sitvfied[3][0]
Bxy_sitv_max = sitvfied[3][1]
Bxy_sitv_mean = sitvfied[3][2]


MVAed = mmcal_functions.MVA([sitv_midpoints_secs,[Bx_sitv,By_sitv,Bz_sitv],[Bxy_sitv_min,Bxy_sitv_max,Bxy_sitv_mean]])
phi_PN16 = MVAed[0][0]
theta_D_PN16 = MVAed[0][1]
theta_B_PN16 = MVAed[0][2]
Bxy_fluct_PN16 = MVAed[0][3]

B_x1_angle = MVAed[1][0]
theta_D = MVAed[1][1]
phi_D = MVAed[1][2]
lam1_lam2 = MVAed[1][3]
lam3_lam2 = MVAed[1][4]

O_z_unfiltered = MVAed[2]
    

################################################################################
################################END OF MVA######################
################################################################################


Mmode_count_PN16 = 0
Mmode_count_SLD08 = 0

Mmode_indices_PN16 = []
Mmode_indices_SLD08 = []
Mmode_indices_overlapping = []
Mmode_indices_PN16_only = []
Mmode_indices_SLD08_only = []




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
starttim=[]

for i in range(0,len(sitv_midpoints_secs)):
    PN16_satisfy=False
    SLD08_satisfy=False
    
    if Bxy_fluct_PN16[i]>C_xy:
        if phi_PN16[i] < C_phi:
            if abs(theta_B_PN16[i]) < C_B:
                if abs(theta_D_PN16[i]) < C_D:
                    Mmode_indices_PN16.append(i)
                    Mmode_count_PN16+=1
                    PN16_satisfy=True
                    
    if lam1_lam2[i] > C_lam12:
        if lam3_lam2[i] > C_lam32:
            #if B_x1_angle[i] < C_MV_B:
            Mmode_indices_SLD08.append(i)
            Mmode_count_SLD08+=1
            SLD08_satisfy=True
            
    if PN16_satisfy:
        O_z_overall_PN16.append(O_z_unfiltered[i])
        
        d_thetaD.append(np.arctan( (1/lam1_lam2[i]) ))
        d_b.append(abs(B_magn[i])*d_g+d_bn)
        #d_b.append(abs(Bz_sitv_mean[i])*d_g+d_bn)
        bxymirror.append(Bxy_sitv_mean[i])
        bzmirror.append(Bz_sitv_mean[i])
        theta_mb.append(theta_B_PN16[i])
        theta_md.append(theta_D_PN16[i])
        test.append(Bxy_sitv_mean[i]*(np.tan(theta_B_PN16[i])-np.tan(theta_D_PN16[i])))
        starttim.append(math.floor(subintervals[i][0]/(3600*24)))
         #O_z.append(Bz_sitv_mean[i]-x1[2]/x1_xy*Bxy_sitv_mean[i])
         
        d_thetaB.append((d_b[-1]/(1+(bzmirror[-1]/bxymirror[-1])**2))*np.sqrt((1/bxymirror[-1])**2+(bzmirror[-1]/bxymirror[-1]**2)**2))   
        error_O_z_PN16.append(np.sqrt(((np.tan(theta_mb[-1])-np.tan(theta_md[-1]))*d_b[-1])**2+((bxymirror[-1]*d_thetaB[-1])/(np.cos(theta_mb[-1]))**2)**2+((bxymirror[-1]*d_thetaD[-1])/(np.cos(theta_md[-1]))**2)**2))

         
        if SLD08_satisfy:
            if len(Mmode_indices_PN16)>0 and len(Mmode_indices_SLD08)>0:
                if Mmode_indices_PN16[-1]==Mmode_indices_SLD08[-1]:
                    Mmode_indices_overlapping.append(i)  
                    O_z_overall_overlapping.append(O_z_overall_PN16[-1])
                    error_O_z_overlapping.append(error_O_z_PN16[-1])
        else:
            Mmode_indices_PN16_only.append(i)
    else:
        if SLD08_satisfy:
            Mmode_indices_SLD08_only.append(i)
        #else: nothing is satisfied
            
        
error=[]
for i in range(len(d_b)):
    d_thetaB.append((d_b[i]/(1+(bzmirror[i]/bxymirror[i])**2))*np.sqrt((1/bxymirror[i])**2+(bzmirror[i]/bxymirror[i]**2)**2))   
    error.append(np.sqrt(((np.tan(theta_mb[i])-np.tan(theta_md[i]))*d_b[i])**2+((bxymirror[i]*d_thetaB[i])/(np.cos(theta_mb[i]))**2)**2+((bxymirror[i]*d_thetaD[i])/(np.cos(theta_md[i]))**2)**2))                 
test=np.array(test)

            
error=np.array(error)
plt.figure()
plt.plot(error,O_z_overall_PN16,'x')
x=np.linspace(0,20,1000)
y=error+5
y1=-error+5
#y=error
#y1=-error
plt.plot(error,y,'r')
plt.ylim(-30,30)
plt.plot(error,y1,'r')

width=plt.hist(error,bins='auto')[1].max()
w=[]
for i in range(len(error)):
    w.append(np.exp(-error[i]**2/(2*width**2)))



dailyO=[]
dailystd=[]
N=[]
sigma=np.array([1.]*len(O_z_overall_PN16))
for i in range(len(mmcal_functions.dailyOff(O_z_overall_PN16,starttim)[0])):
   dailyO.append( mmcal_functions.wadapkde(mmcal_functions.dailyOff(O_z_overall_PN16,starttim)[0][i],w))
   dailystd.append(np.std(mmcal_functions.dailyOff(O_z_overall_PN16,starttim)[0][i]))
   N.append(len(mmcal_functions.dailyOff(O_z_overall_PN16,starttim)[0][i]))
   
dailystd=np.array(dailystd)   
N=np.array(N)





################################################################################
    
######################## PLOTTING ##############################################

################################################################################
Mmode_count_overlapping = len(Mmode_indices_overlapping)
sitv_midpoints_days = np.array(sitv_midpoints_secs)/3600/24

#sitv_midpoints_days = (subintervals[:,0]+subintervals[:,1])/2/3600/24

Mmode_times = mmcal_functions.get_Mmode_times(sitv_midpoints_days,Mmode_indices_PN16_only,Mmode_indices_SLD08_only,Mmode_indices_overlapping)
Mmode_sitv_times_PN16_only = Mmode_times[0]
Mmode_sitv_times_SLD08_only = Mmode_times[1]
Mmode_sitv_times_overlapping = Mmode_times[2]


print("PN16:",Mmode_count_PN16)
print("SLD08:",Mmode_count_SLD08)
print("overlapping:",len(Mmode_indices_overlapping))
print("PN16 only:",len(Mmode_sitv_times_PN16_only))
print("SLD08 only:",len(Mmode_sitv_times_SLD08_only))    




#Obtaining the pdf using kde
O_z_PN16_T = np.transpose(np.array([O_z_overall_PN16]))
kde= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(O_z_PN16_T)
log_pdf = kde.score_samples(O_z_PN16_T)
pdf = np.exp(log_pdf)

        
O_z_PN16_only = []
O_z_overlapping = []
pdf_PN16_only = []
pdf_overlapping1 = []

i=0
j=0
k=0
#while i < len(subintervals):
while j+k < len(O_z_overall_PN16):
    if j<len(Mmode_indices_PN16_only):
        if i==Mmode_indices_PN16_only[j]:
            O_z_PN16_only.append(O_z_overall_PN16[j+k])
            pdf_PN16_only.append(pdf[j+k])
            j+=1
    if k<len(Mmode_indices_overlapping):
        if i==Mmode_indices_overlapping[k]:
            O_z_overlapping.append(O_z_overall_PN16[j+k])
            pdf_overlapping1.append(pdf[j+k])
            k+=1
    i+=1


kde_overlapping2= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(np.transpose(np.array([O_z_overlapping])))
log_pdf_overlapping2 = kde_overlapping2.score_samples(np.transpose(np.array([O_z_overlapping])))
pdf_overlapping2 = np.exp(log_pdf_overlapping2)

################################################################################
    
######################## PRINTING ##############################################

################################################################################

print("O_z PN16:",O_z_overall_PN16[np.argmax(pdf)],np.mean(O_z_overall_PN16),np.std(O_z_overall_PN16))
print("O_z both:",O_z_overlapping[np.argmax(pdf_overlapping2)],np.mean(O_z_overlapping),np.std(O_z_overlapping))



################################################################################
    
######################## PLOTTING ##############################################

################################################################################



error_O_z_PN16 = np.array(error_O_z_PN16)
error_O_z_overlapping = np.array(error_O_z_overlapping)

width_PN16=plt.hist(error_O_z_PN16,bins='auto')[1].max()
w_PN16=[]
for i in range(len(error_O_z_PN16)):
    w_PN16.append(np.exp(-error_O_z_PN16[i]**2/(2*width_PN16**2)))  
    
width_overlapping=plt.hist(error_O_z_overlapping,bins='auto')[1].max()
w_overlapping=[]
for i in range(len(error_O_z_overlapping)):
    w_overlapping.append(np.exp(-error_O_z_overlapping[i]**2/(2*width_overlapping**2)))  


h_PN16=bandwidth.sskernel(np.array(O_z_overall_PN16))[2]
x_PN16=np.linspace(-20,20,100000)
ker=0
for i in range(len(O_z_overall_PN16)):
    mu=O_z_overall_PN16
    y=1/(h_PN16 * np.sqrt(2 * np.pi)*w_PN16[i])*np.exp( - (x_PN16 - mu[i])**2 / (2 * h_PN16**2))
    ker+=y
ker_PN16=ker/len(O_z_overall_PN16)

h_overlapping=bandwidth.sskernel(np.array(O_z_overall_overlapping))[2]
x_overlapping=np.linspace(-20,20,100000)
ker=0
for i in range(len(O_z_overall_overlapping)):
    mu=O_z_overall_overlapping
    y=1/(h_overlapping * np.sqrt(2 * np.pi)*w_overlapping[i])*np.exp( - (x_overlapping - mu[i])**2 / (2 * h_overlapping**2))
    ker+=y
ker_overlapping=ker/len(O_z_overall_overlapping)


kde_overall_PN16= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(np.transpose(np.array([O_z_overall_PN16])))
log_pdf_overall_PN16 = kde_overall_PN16.score_samples(np.transpose(np.array([O_z_overall_PN16])))
pdf_overall_PN16 = np.exp(log_pdf_overall_PN16)

kde_overall_overlapping= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(np.transpose(np.array([O_z_overall_overlapping])))
log_pdf_overall_overlapping = kde_overall_overlapping.score_samples(np.transpose(np.array([O_z_overall_overlapping])))
pdf_overall_overlapping = np.exp(log_pdf_overall_overlapping)

print("FINAL SUMMARY:\n")
print("PN16:",len(O_z_overall_PN16))
print("SLD08:",len(Mmode_indices_SLD08))
print("overlapping:",len(O_z_overall_overlapping))
print("PN16 only:",len(Mmode_indices_PN16_only))
print("SLD08 only:",len(Mmode_indices_SLD08_only))

print("O_z PN16:",O_z_overall_PN16[np.argmax(pdf_overall_PN16)],np.mean(O_z_overall_PN16),np.std(O_z_overall_PN16))
print("O_z PN16 wadapkde:",x_PN16[np.where(ker_PN16==np.max(ker_PN16))[0][0]])
print("O_z both:",O_z_overall_overlapping[np.argmax(pdf_overall_overlapping)],np.mean(O_z_overall_overlapping),np.std(O_z_overall_overlapping))
print("O_z both wadapkde:",x_overlapping[np.where(ker_overlapping==np.max(ker_overlapping))[0][0]])



if PLOT:
    


    f1=plt.figure()
    ax1 = f1.add_subplot(111)
    
    ax1.plot_date(t_days,B_mag,fmt='-',linewidth=1.0,color='black')
    ax1.set_title(r"B-field magnitude time series $(t_{si}=%d,t_{sh}=%d)$"%(t_int,shift))
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
    ax1.axvline(Mmode_sitv_times_PN16_only[0],alpha=0.4,color='orange',linewidth=1.0,label='PN16 only')
    if len(Mmode_sitv_times_SLD08_only) > 0:
        ax1.axvline(Mmode_sitv_times_SLD08_only[0],alpha=0.4,color='green',linewidth=1.0,label='SLD08(no[4]) only')    
    ax1.axvline(Mmode_sitv_times_overlapping[0],alpha=0.4,color='blue',linewidth=1.0,label='PN16&SLD08(no[4])')    
    ax1.legend()

    """
    f2=plt.figure()
    ax21 = f2.add_subplot(311)
    ax22 = f2.add_subplot(312)
    ax23 = f2.add_subplot(313)
    
    ax21.plot(t_days,B_x,linewidth=1.0)
    ax21.plot_date(sitv_midpoints_days,Bx_sitv_mean,fmt='-',linewidth=1.0)
    ax21.set_title("Cartesian B-field time series")
    ax21.set_ylabel(r"$B_{x}$ (nT)")
    ax22.plot(t_days,B_y,linewidth=1.0)
    ax22.plot_date(sitv_midpoints_days,By_sitv_mean,fmt='-',linewidth=1.0)
    ax22.set_ylabel(r"$B_{y}$ (nT)")
    ax23.plot(t_days,B_z,linewidth=1.0)
    ax23.plot_date(sitv_midpoints_days,Bz_sitv_mean,fmt='-',linewidth=1.0)
    ax23.set_ylabel(r"$B_{z}$ (nT)")
    ax23.set_xlabel("Time")
    """

    
    f3=plt.figure()
    ax31 = f3.add_subplot(311)
    ax32 = f3.add_subplot(312)
    ax33 = f3.add_subplot(313)
    
    ax31.plot_date(t_days,B_x/B_mag,fmt='-',linewidth=1.0)
    ax31.set_ylabel(r"$B_{x}$ (nT)")
    ax31.set_title("Normalised Cartesian B-field time series")
    ax32.plot_date(t_days,B_y/B_mag,fmt='-',linewidth=1.0)
    ax32.set_ylabel(r"$B_{y}$ (nT)")
    ax33.plot_date(t_days,B_z/B_mag,fmt='-',linewidth=1.0)
    ax33.set_ylabel(r"$B_{z}$ (nT)")
    ax33.set_xlabel("Time")

    """
    f4=plt.figure()
    ax40 = f4.add_subplot(311)
    ax41 = f4.add_subplot(312)
    ax42 = f4.add_subplot(313)    
    
    ax40.plot_date(t_days,B_mag,fmt='-',linewidth=1.0)
    ax41.plot_date(t_days,np.array(theta)*180/np.pi,fmt='-',linewidth=1.0)
    ax42.plot_date(t_days,np.array(phi)*180/np.pi,fmt='-',linewidth=1.0)
    ax40.set_title("Sph. Polar B-field time series")
    ax40.set_ylabel(r"$B_{mag}$ (nT)")
    ax41.set_ylabel(r"$\theta$ (degs)")
    ax42.set_ylabel(r"$\phi$ (degs)")
    ax42.set_xlabel("Time")
    """

    """
    f5=plt.figure()
    ax5 = f5.add_subplot(111)
    
    ax5.plot_date(sitv_midpoints_days,np.array(B_x1_angle)*180/np.pi,fmt='-',linewidth=1.0)
    ax5.plot_date(sitv_midpoints_days,np.array(len(sitv_midpoints_days)*[30]),fmt='-',linewidth=1.0)
    ax5.set_ylabel("B_D_angle (degs)")
    ax5.set_title("Time Series of Angle between MV dir. and B-field dir.")
    ax5.set_xlabel("Time")
    """
    
    f6=plt.figure()
    ax61 = f6.add_subplot(211)
    ax62 = f6.add_subplot(212)
    
    ax61.plot_date(sitv_midpoints_days,np.abs(np.array(phi_PN16))*180/np.pi,fmt='-',linewidth=1.0)
    ax62.plot_date(sitv_midpoints_days,np.array(theta_D_PN16)*180/np.pi,fmt='-',linewidth=1.0)
    ax61.plot_date(sitv_midpoints_days,[20]*len(phi_PN16),fmt='-',linewidth=1.0)
    ax62.plot_date(sitv_midpoints_days,[30]*len(theta_D_PN16),fmt='-',linewidth=1.0)
    ax62.plot_date(sitv_midpoints_days,[-30]*len(theta_D_PN16),fmt='-',linewidth=1.0)
    ax61.set_ylabel(r"$\phi$")
    ax62.set_ylabel(r"$\theta_{D}$")
    
    f7=plt.figure()
    ax71 = f7.add_subplot(211)
    ax72 = f7.add_subplot(212)   
    
    ax71.plot_date(sitv_midpoints_days,np.array(theta_B_PN16)*180/np.pi,fmt='-',linewidth=1.0)
    ax71.plot_date(sitv_midpoints_days,[30]*len(theta_B_PN16),fmt='-',linewidth=1.0)
    ax71.plot_date(sitv_midpoints_days,[-30]*len(theta_B_PN16),fmt='-',linewidth=1.0)    
    ax72.plot_date(sitv_midpoints_days,Bxy_fluct_PN16,fmt='-',linewidth=1.0)    
    ax72.plot_date(sitv_midpoints_days,[C_xy]*len(Bxy_fluct_PN16),fmt='-',linewidth=1.0)          
    ax71.set_ylabel(r"$\theta_{B}$")
    ax72.set_ylabel(r"$\frac{\delta B_{xy}}{B_{xy}}$")
    
    """
    f8=plt.figure()
    ax81 = f8.add_subplot(211)
    ax82 = f8.add_subplot(212)          
    ax81.plot_date(sitv_midpoints_days,np.array(lam1_lam2),fmt='-',linewidth=1.0)
    ax81.plot_date(sitv_midpoints_days,[C_lam12]*len(lam1_lam2),fmt='-',linewidth=1.0) 
    ax82.plot_date(sitv_midpoints_days,np.array(lam3_lam2),fmt='-',linewidth=1.0)    
    ax82.plot_date(sitv_midpoints_days,[C_lam32]*len(lam3_lam2),fmt='-',linewidth=1.0)   
    ax81.set_ylabel(r"$\frac{\lambda_{1}}{\lambda_{2}}$")
    ax82.set_ylabel(r"$\frac{\lambda_{3}}{\lambda_{2}}$")
    """
    
    f9=plt.figure()
    ax9 = f9.add_subplot(111)
    ax9.scatter(O_z_overall_PN16,pdf_overall_PN16,s=1.0,color='orange',label="PN16")
    ax9.plot(x_PN16,ker_PN16,linewidth=1.0,color='red',label='adaptive weighted kde')
    ax9.set_title("PDF from KDE: C3 "+plot_time_text+", %d PN16 intervals, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(len(O_z_overall_PN16),t_int,shift,artificial_Bz_offset))
    ax9.set_ylabel("Prob. Density")
    ax9.set_xlabel(r"$O_{z}$ (nT)")
    ax9.legend()        

    f10=plt.figure()
    ax10 = f10.add_subplot(111)
    ax10.scatter(O_z_overall_overlapping,pdf_overall_overlapping,s=1.0,color='blue',label="PN and SLD(2)")
    ax10.plot(x_overlapping,ker_overlapping,linewidth=1.0,color='red',label='adaptive weighted kde')
    ax10.set_title("PDF from KDE: C3 "+plot_time_text+", %d PN&SLD intervals, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(len(O_z_overall_overlapping),t_int,shift,artificial_Bz_offset))
    ax10.set_ylabel("Prob. Density")
    ax10.set_xlabel(r"$O_{z}$ (nT)")
    ax10.legend()     
    

    plt.show()        






