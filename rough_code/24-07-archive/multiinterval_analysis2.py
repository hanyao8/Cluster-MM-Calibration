# -*- coding: utf-8 -*-
"""
Created on Sat Jul  7 18:11:14 2018

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
import math
import bandwidth
import mmcal_functions

#import pywt
#from scipy.integrate import quad

from sklearn.neighbors.kde import KernelDensity

#offset_analysis
#contains code from 'mmode_search' part of main.py
#dsl based data must be used for this to give valid results

################################################################################

######################## UI AND SETTINGS ##############################################

################################################################################

ARTIFICIAL_OFFSET = True
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
csv_file_name = "C3_CP_FGM_5VPS__200511_sheathselected"

sheath_file_name = "C3_CP_sheathintervals__200511"

#data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T00:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T22:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T21:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T22:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T00:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T23:59:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))

yyyy_s = '2005'
mm_s = '11'
dd_s = '01'
hh24_s = '00'

yyyy_e = '2005'
mm_e = '11'
dd_e = '30'
hh24_e = '23'

#plot_time_text = yyyy_s+'/'+mm_s+'/'+dd_s+' '+hh24_s+'~'+hh24_e
plot_time_text = yyyy_s+'/'+mm_s+'/'+dd_s+' '+hh24_s+'~'+ dd_e+' ' +hh24_e

overall_data_start_time = matplotlib.dates.date2num(datetime.strptime(yyyy_s+'-'+mm_s+'-'+dd_s+'T'+hh24_s+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
overall_data_end_time = matplotlib.dates.date2num(datetime.strptime(yyyy_e+'-'+mm_e+'-'+dd_e+'T'+hh24_e+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))

#plot_time_text = "Jul 2008"
plot_time_text = "Nov 2005"

t_int = 180
shift = 10

C_xy = 0.3
C_B = 30*np.pi/180
C_D = 30*np.pi/180
C_phi = 20*np.pi/180

C_lam12 = 1.5
C_lam32 = 0.3
#C_MV_B = 30*np.pi/180

var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']
#var_names = ['Time','Half Interval','Bx','By','Bz','Bt']

#csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ csv_file_name + ".csv",names=var_names)
csv_df = pd.read_csv(os.getcwd()+"\\"+ csv_file_name + ".csv",names=var_names)
csv_df.head()
df_arr = csv_df.values
print("Progress: Main DSL data file loaded")

#csv_df = pd.read_csv(os.getcwd()+"\\data\\"+ sheath_file_name + ".csv")
sh_var_names = ['sheath start','sheath end']
csv_df = pd.read_csv(os.getcwd()+"\\"+ sheath_file_name + ".csv",names=sh_var_names)
csv_df.head()
sheath_df_arr = csv_df.values
print("Progress: Sheath-pass file loaded")

#sheath_df_arr = sheath_df_arr[1:2]

################################################################################

######################## FUNCTIONS #############################################

################################################################################

if csv_file_name[:3]!='DSL':
    print("Ensure B-field data is in DSL coordinates")

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


t_overall = df_arr[:,0]
B_x_overall = df_arr[:,2]
B_y_overall = df_arr[:,3]
B_z_overall = df_arr[:,4]
rangelist_overall = (df_arr[:,9]).astype(int)
print(len(B_z_overall))


    

if simulate_spinres:
    data_overall_spinres = fivevps_2_spinres(t_overall,B_x_overall,B_y_overall,B_z_overall)
    t_overall = data_overall_spinres[0]
    B_x_overall = np.array(data_overall_spinres[1])
    B_y_overall = np.array(data_overall_spinres[2])
    B_z_overall = np.array(data_overall_spinres[3])          


#B_mag = df_arr[:,5]
print("Progress: Variables extracted from dataframe array")

#!!!!!!!!!!!!! artificial offset set here:
if ARTIFICIAL_OFFSET:
    B_z_overall = np.array(len(B_z_overall)*[artificial_Bz_offset]) + np.array(B_z_overall)

"""
if SETTING=='themis_validation':
    t_days = t

else:
"""    
""" 
t_datetime = []
for i in range(0,len(t_overall)):
    strpdtime=datetime.strptime(t_overall[i],'%Y-%m-%dT%H:%M:%S.%fZ')
    t_datetime.append(strpdtime)
t_days_overall = matplotlib.dates.date2num(t_datetime)
"""
t_days_overall=[]
for i in range(0,len(t_overall)):
    t_days_overall.append(matplotlib.dates.date2num(datetime.strptime(t_overall[i],'%Y-%m-%dT%H:%M:%S.%fZ')))
t_days_overall=np.array(t_days_overall)
print("Progress: Time parsing finished, big loop about to begin")


#############################################################
################### BIG LOOP ################################
#############################################################

count_stats = [[],[],[],[],[]]
kde_stats = [[[],[],[]],[[],[],[]]]
sheath_itv_midpoints = []
sheath_itv_midpoints_plot1 = []
O_z_overall_PN16 = []
O_z_overall_overlapping = []
error_O_z_PN16 = []
error_O_z_overlapping = []

range_sitv_PN16 = []
O_z_overall_PN16_r2 = []
O_z_overall_PN16_r3 = []
O_z_overall_PN16_r23 = []


for sheathpass_num in range(0,len(sheath_df_arr)):
    print("Progress: Sheath-pass:",sheathpass_num)
    
    ts_str = sheath_df_arr[sheathpass_num][0]
    tf_str = sheath_df_arr[sheathpass_num][1]
    
    ts_num = matplotlib.dates.date2num(datetime.strptime(ts_str,'%Y-%m-%dT%H:%M:%S.%fZ'))
    tf_num = matplotlib.dates.date2num(datetime.strptime(tf_str,'%Y-%m-%dT%H:%M:%S.%fZ'))
    sheath_itv_midpoints.append((ts_num+tf_num)/2)

    if ts_num > t_days_overall[0]:
        data_start_index = np.argmax(t_days_overall>ts_num)
    else:
        data_start_index = 0
    if tf_num < t_days_overall[-1]:
        data_end_index = np.argmax(t_days_overall>tf_num)
    else:
        data_end_index = len(t_days_overall)
    
    print("Progress: 1")
    
    t = t_overall[data_start_index:data_end_index]
    B_x = B_x_overall[data_start_index:data_end_index]
    B_y = B_y_overall[data_start_index:data_end_index]
    B_z = B_z_overall[data_start_index:data_end_index]
    rangelist = rangelist_overall[data_start_index:data_end_index]
    #B_mag = B_mag[data_start_index:data_end_index]
    B_x.astype(float)
    B_y.astype(float)
    B_z.astype(float)
    B_mag = (B_x**2+B_y**2+B_z**2)**(0.5)
    

    
    t_days = t_days_overall[data_start_index:data_end_index]
    t_secs = t_days*3600*24

    print("Progress: 2")    
    
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
        
    print("Progress: 3")    
    
    sitvfied = mmcal_functions.sitvfy([t_secs,[B_x,B_y,B_z],B_mag,rangelist],[t_int,shift])
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
    range_sitv = sitvfied[4]
    B_magn=[np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2) for i in range(0,len(Bx_sitv_mean)) ]

    print("Progress: 4")

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

    print("Progress: 5")
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
    
    
    O_z_PN16 = []
    
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
    #using PN16 criteria 
    starttim=[]
    error=[]
    
    for i in range(0,len(sitv_midpoints_secs)):

        PN16_satisfy=False
        SLD08_satisfy=False
        

        #if Bmag_sitv_mean > 10:
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

            O_z_PN16.append(O_z_unfiltered[i])
            O_z_overall_PN16.append(O_z_unfiltered[i])

            d_thetaD.append(np.arctan(1/lam1_lam2[i]))
            d_b.append(abs(B_magn[i])*d_g+d_bn)
            #d_b.append(abs(Bz_sitv_mean[i])*d_g+d_bn)
            bxymirror.append(Bxy_sitv_mean[i])
            bzmirror.append(Bz_sitv_mean[i])
            theta_mb.append(theta_B_PN16[i])
            theta_md.append(theta_D_PN16[i])
            test.append(Bxy_sitv_mean[i]*(np.tan(theta_B_PN16[i])-np.tan(theta_D_PN16[i])))
            starttim.append(math.floor(subintervals[i][0]/(60**2*24)))

            d_thetaB.append((d_b[-1]/(1+(bzmirror[-1]/bxymirror[-1])**2))*np.sqrt((1/bxymirror[-1])**2+(bzmirror[-1]/bxymirror[-1]**2)**2))   
            error.append(np.sqrt(((np.tan(theta_mb[-1])-np.tan(theta_md[-1]))*d_b[-1])**2+((bxymirror[-1]*d_thetaB[-1])/(np.cos(theta_mb[-1]))**2)**2+((bxymirror[-1]*d_thetaD[-1])/(np.cos(theta_md[-1]))**2)**2))         
            error_O_z_PN16.append(error[-1])
            
            range_sitv_PN16.append(range_sitv[i])
            if range_sitv[i]==2:
                O_z_overall_PN16_r2.append(O_z_PN16[-1])
            elif range_sitv[i]==3:
                O_z_overall_PN16_r3.append(O_z_PN16[-1])      
            elif range_sitv[i]==23:
                O_z_overall_PN16_r23.append(O_z_PN16[-1])      
            

        if PN16_satisfy and SLD08_satisfy:
            if len(Mmode_indices_PN16)>0 and len(Mmode_indices_SLD08)>0:
                if Mmode_indices_PN16[-1]==Mmode_indices_SLD08[-1]:
                    Mmode_indices_overlapping.append(i)    
                    O_z_overall_overlapping.append(O_z_unfiltered[i])
                    error_O_z_overlapping.append(error[-1])
            
        if PN16_satisfy and not(SLD08_satisfy):
            Mmode_indices_PN16_only.append(i)
        if SLD08_satisfy and not(PN16_satisfy):
            Mmode_indices_SLD08_only.append(i)
    

    #test=np.array(test)    
    #error=np.array(error)
    
    print("Progress: 6")  
    ################################################################################
    
    ########################          ##############################################
    
    ################################################################################
    Mmode_count_overlapping = len(Mmode_indices_overlapping)
    
    O_z_PN16 = np.array(O_z_PN16)
    
    
    
    
    sitv_midpoints_days = (subintervals[:,0]+subintervals[:,1])/2/3600/24
    Mmode_sitv_times_PN16 = []
    Mmode_sitv_times_SLD08 = []
    Mmode_sitv_times_PN16_days = []
    Mmode_sitv_times_SLD08_days = []
    Mmode_sitv_times_PN16_only = []
    Mmode_sitv_times_SLD08_only = []
    
    for i in range(0,len(Mmode_indices_PN16)):
        Mmode_sitv_times_PN16.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_PN16[i]]))
        Mmode_sitv_times_PN16_days.append(sitv_midpoints_days[Mmode_indices_PN16[i]])
    for i in range(0,len(Mmode_indices_SLD08)):
        Mmode_sitv_times_SLD08.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_SLD08[i]]))   
        Mmode_sitv_times_SLD08_days.append(sitv_midpoints_days[Mmode_indices_SLD08[i]])
    for i in range(0,len(Mmode_indices_PN16_only)):
        Mmode_sitv_times_PN16_only.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_PN16_only[i]]))
    for i in range(0,len(Mmode_indices_SLD08_only)):
        Mmode_sitv_times_SLD08_only.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_SLD08_only[i]]))   
    
    
    Mmode_sitv_times_overlapping = [] 
    Mmode_sitv_times_overlapping_days = [] 
    for i in range(0,len(Mmode_indices_overlapping)):        
        Mmode_sitv_times_overlapping.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_overlapping[i]]))
        Mmode_sitv_times_overlapping_days.append(sitv_midpoints_days[Mmode_indices_overlapping[i]])
    
    print("PN16:",Mmode_count_PN16)
    print("SLD08:",Mmode_count_SLD08)
    print("overlapping:",len(Mmode_indices_overlapping))
    print("PN16 only:",len(Mmode_sitv_times_PN16_only))
    print("SLD08 only:",len(Mmode_sitv_times_SLD08_only))

        
    if sitv_in_sheath_frac_analysis:
        sheath_start_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T06:58:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
        sheath_end_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T11:40:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
        
        PN16_frac_in_sheath = (np.argmax(np.array(Mmode_sitv_times_PN16_days)>sheath_end_time) - np.argmax(np.array(Mmode_sitv_times_PN16_days)>sheath_start_time) )/Mmode_count_PN16
        SLD08_frac_in_sheath = (np.argmax(np.array(Mmode_sitv_times_SLD08_days)>sheath_end_time) - np.argmax(np.array(Mmode_sitv_times_SLD08_days)>sheath_start_time) )/Mmode_count_SLD08  
        overlapping_frac_in_sheath = (np.argmax(np.array(Mmode_sitv_times_overlapping_days)>sheath_end_time) - np.argmax(np.array(Mmode_sitv_times_overlapping_days)>sheath_start_time) )/len(Mmode_indices_overlapping)
        
    
        
        print(PN16_frac_in_sheath)
        print(SLD08_frac_in_sheath)
        print(overlapping_frac_in_sheath)
    
    
    #if len(O_z) > 0:
    if len(Mmode_indices_overlapping) > 0:
        
        sheath_itv_midpoints_plot1.append((ts_num+tf_num)/2)
        
        #Obtaining the pdf using kde
        O_z_T = np.transpose(np.array([O_z_PN16]))
        kde= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(O_z_T)
        log_pdf = kde.score_samples(O_z_T)
        pdf = np.exp(log_pdf)
        
                
        O_z_PN16_only = []
        O_z_overlapping = []
        pdf_PN16_only = []
        pdf_overlapping1 = []
        
        i=0
        j=0
        k=0
        #while i < len(subintervals):
        while j+k < len(O_z_PN16):
            if j<len(Mmode_indices_PN16_only):
                if i==Mmode_indices_PN16_only[j]:
                    O_z_PN16_only.append(O_z_PN16[j+k])
                    pdf_PN16_only.append(pdf[j+k])
                    j+=1
            if k<len(Mmode_indices_overlapping):
                if i==Mmode_indices_overlapping[k]:
                    O_z_overlapping.append(O_z_PN16[j+k])
                    pdf_overlapping1.append(pdf[j+k])
                    k+=1
            i+=1
        """
        for i in range(0,len(Mmode_indices_PN16_only)):
            O_z_PN16_only.append(O_z_wfill[Mmode_indices_PN16_only[i]] )
            pdf_PN16_only.append(pdf_wfill[Mmode_indices_PN16_only[i]] )
        for i in range(0,len(Mmode_indices_overlapping)):
            O_z_overlapping.append(O_z_wfill[Mmode_indices_overlapping[i]] )
            pdf_overlapping.append(pdf_wfill[Mmode_indices_overlapping[i]] )
        """
        
        
        kde_overlapping2= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(np.transpose(np.array([O_z_overlapping])))
        log_pdf_overlapping2 = kde_overlapping2.score_samples(np.transpose(np.array([O_z_overlapping])))
        pdf_overlapping2 = np.exp(log_pdf_overlapping2)
        
        ################################################################################
            
        ######################## PRINTING ##############################################
        
        ################################################################################
        
        print("O_z PN16:",O_z_PN16[np.argmax(pdf)],np.mean(O_z_PN16),np.std(O_z_PN16))
        print("O_z both:",O_z_overlapping[np.argmax(pdf_overlapping2)],np.mean(O_z_overlapping),np.std(O_z_overlapping))
        kde_stats[0][0].append(O_z_PN16[np.argmax(pdf)])
        kde_stats[0][1].append(np.mean(O_z_PN16))
        kde_stats[0][2].append(np.std(O_z_PN16))
        kde_stats[1][0].append(O_z_overlapping[np.argmax(pdf_overlapping2)])
        kde_stats[1][1].append(np.mean(O_z_overlapping))
        kde_stats[1][2].append(np.std(O_z_overlapping))
            
        count_stats[0].append(Mmode_count_PN16)
        count_stats[1].append(Mmode_count_SLD08)
        count_stats[2].append(len(Mmode_indices_overlapping))
        count_stats[3].append(len(Mmode_sitv_times_PN16_only))
        count_stats[4].append(len(Mmode_sitv_times_SLD08_only))        

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

count_stats=np.array(count_stats)
kde_stats=np.array(kde_stats)

kde_overall_PN16= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(np.transpose(np.array([O_z_overall_PN16])))
log_pdf_overall_PN16 = kde_overall_PN16.score_samples(np.transpose(np.array([O_z_overall_PN16])))
pdf_overall_PN16 = np.exp(log_pdf_overall_PN16)

kde_overall_overlapping= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(np.transpose(np.array([O_z_overall_overlapping])))
log_pdf_overall_overlapping = kde_overall_overlapping.score_samples(np.transpose(np.array([O_z_overall_overlapping])))
pdf_overall_overlapping = np.exp(log_pdf_overall_overlapping)

print("\nFINAL SUMMARY:\n")
print("PN16:",sum(count_stats[0]))
print("SLD08:",sum(count_stats[1]))
print("overlapping:",sum(count_stats[2]))
print("PN16 only:",sum(count_stats[3]))
print("SLD08 only:",sum(count_stats[4]))

print("$O_{z}$ PN16:")
print("$O_{zf}=$", O_z_overall_PN16[np.argmax(pdf_overall_PN16)] )
print("$\langle O_{z} \rangle=$", np.mean(O_z_overall_PN16) )
print("$\sigma_{Oz}=$", np.std(O_z_overall_PN16) ) 
print("$O_{z}$ PN16 wadapkde:")
print("$O_{zf}=$" ,x_PN16[np.where(ker_PN16==np.max(ker_PN16))[0][0]] )

print("$O_{z}$ both:")
print("$O_{zf}=$", O_z_overall_overlapping[np.argmax(pdf_overall_overlapping)] )
print("$\langle O_{z} \rangle=$", np.mean(O_z_overall_overlapping) )
print("$\sigma_{Oz}=$", np.std(O_z_overall_overlapping) )
print("$O_{z}$ both wadapkde:")
print("$O_{zf}=$" ,x_overlapping[np.where(ker_overlapping==np.max(ker_overlapping))[0][0]]) 

if PLOT:
    
    sitv_midpoints = (subintervals[:,0]+subintervals[:,1])/2/3600/24
    
    f1=plt.figure()
    ax11 = f1.add_subplot(211)
    ax12 = f1.add_subplot(212)    
    ax11.plot_date(sheath_itv_midpoints_plot1,kde_stats[0][0],markersize=2.0)
    for i, anno in enumerate(count_stats[0]):
        ax11.annotate(anno,(sheath_itv_midpoints_plot1[i],kde_stats[0][0][i]))
    ax12.plot_date(sheath_itv_midpoints_plot1,kde_stats[0][1],markersize=2.0)
    for i, anno in enumerate(count_stats[0]):
        ax12.annotate(anno,(sheath_itv_midpoints_plot1[i],kde_stats[0][1][i]))
    ax11.set_title(r"$O_{zf}$ estimates per sheath-pass C3 "+plot_time_text+", PN16, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(t_int,shift,artificial_Bz_offset))
    #ax11.set_title(r"$O_{zf}$ estimates per day C3 "+plot_time_text+", PN16, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(t_int,shift,artificial_Bz_offset))    
    ax11.set_ylabel("KDE dist maximum")
    ax12.set_ylabel("KDE dist mean")
    
    f2=plt.figure()
    ax21 = f2.add_subplot(211)
    ax22 = f2.add_subplot(212)
    ax21.plot_date(sheath_itv_midpoints_plot1,kde_stats[1][0],markersize=2.0)
    for i, anno in enumerate(count_stats[2]):
        ax21.annotate(anno,(sheath_itv_midpoints_plot1[i],kde_stats[1][0][i]))
    ax22.plot_date(sheath_itv_midpoints_plot1,kde_stats[1][1],markersize=2.0)    
    for i, anno in enumerate(count_stats[2]):
        ax22.annotate(anno,(sheath_itv_midpoints_plot1[i],kde_stats[1][1][i]))
    ax21.set_title(r"$O_{zf}$ estimates per sheath-pass C3 "+plot_time_text+", PN16&2eigen-conditions, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(t_int,shift,artificial_Bz_offset))
    #ax21.set_title(r"$O_{zf}$ estimates per day C3 "+plot_time_text+", PN16&2eigen-conditions, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(t_int,shift,artificial_Bz_offset))    
    ax21.set_ylabel("KDE dist maximum")
    ax22.set_ylabel("KDE dist mean")


    f3=plt.figure()
    ax3 = f3.add_subplot(111)
    ax3.scatter(O_z_overall_PN16,pdf_overall_PN16,s=1.0,color='orange',label="PN16")
    ax3.plot(x_PN16,ker_PN16,linewidth=1.0,color='red',label='adaptive weighted kde')
    ax3.set_title("PDF from KDE: C3 "+plot_time_text+", %d PN16 intervals, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(len(O_z_overall_PN16),t_int,shift,artificial_Bz_offset))
    ax3.set_ylabel("Prob. Density")
    ax3.set_xlabel(r"$O_{z}$ (nT)")
    ax3.legend()

    f4=plt.figure()
    ax4 = f4.add_subplot(111)
    ax4.scatter(O_z_overall_overlapping,pdf_overall_overlapping,s=1.0,color='blue',label="PN and SLD(2)")
    ax4.plot(x_overlapping,ker_overlapping,linewidth=1.0,color='red',label='adaptive weighted kde')
    ax4.set_title("PDF from KDE: C3 "+plot_time_text+", %d PN&SLD intervals, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(len(O_z_overall_overlapping),t_int,shift,artificial_Bz_offset))
    ax4.set_ylabel("Prob. Density")
    ax4.set_xlabel(r"$O_{z}$ (nT)")
    ax4.legend()
    
    f5 = plt.figure()
    ax5 = f5.add_subplot(111)
    ax5.plot(error_O_z_PN16,O_z_overall_PN16,'x')
    #x=np.linspace(0,20,1000)
    ax5.plot(error_O_z_PN16,error_O_z_PN16+5,'r')
    ax5.set_ylim(-30,30)
    ax5.plot(error_O_z_PN16,-error_O_z_PN16+5,'r')
    
    f6 = plt.figure()
    ax6 = f6.add_subplot(111)
    ax6.plot(error_O_z_overlapping,O_z_overall_overlapping,'x')
    #x=np.linspace(0,20,1000)
    ax6.plot(error_O_z_overlapping,error_O_z_overlapping+5,'r')
    ax6.set_ylim(-30,30)
    ax6.plot(error_O_z_overlapping,-error_O_z_overlapping+5,'r')
    
    plt.show()        






