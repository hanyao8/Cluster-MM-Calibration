# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 13:36:49 2018

@author: Choon
"""

#multiinterval_analysis

#import io
import os
import numpy as np
#import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#import time
from datetime import datetime
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

artificial_Bz_offset = 5.0 #nT (in DSL coordinates- hence in spin axis direction)
if not(ARTIFICIAL_OFFSET):
    artificial_Bz_offset = 0.0

#csv_file_name = "DSL_C3_CP_FGM_5VPS__20061103_000000_20061104_000000_V140305"
#csv_file_name = "DSL_C3_CP_FGM_5VPS__20060301_000000_20060302_000000_V140305"
#csv_file_name = "DSL_C3_CP_FGM_5VPS__20060302_000000_20060307_000000_V140305"
#csv_file_name = "DSL_C3_CP_FGM_5VPS__20060307_000000_20060312_000000_V140305"
csv_file_name = "DSL_THEMIS_C_FGM_SPINRES_2008Jul"

sheath_file_name = "THEMIS_C_sheathintervals__2008Jul"

#data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T00:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-03-01T22:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T21:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T22:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_start_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T00:00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
#data_end_time = matplotlib.dates.date2num(datetime.strptime('2006-11-03T23:59:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
"""
yyyy_s = '2006'
mm_s = '03'
dd_s = '01'
hh24_s = '00'

yyyy_e = '2006'
mm_e = '03'
dd_e = '31'
hh24_e = '23'

#plot_time_text = yyyy_s+'/'+mm_s+'/'+dd_s+' '+hh24_s+'~'+hh24_e
plot_time_text = yyyy_s+'/'+mm_s+'/'+dd_s+' '+hh24_s+'~'+ dd_e+' ' +hh24_e

overall_data_start_time = matplotlib.dates.date2num(datetime.strptime(yyyy_s+'-'+mm_s+'-'+dd_s+'T'+hh24_s+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
overall_data_end_time = matplotlib.dates.date2num(datetime.strptime(yyyy_e+'-'+mm_e+'-'+dd_e+'T'+hh24_e+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
"""
plot_time_text = "Jul 2008"

t_int = 180
shift = 10

C_xy = 0.3
C_B = 30*np.pi/180
C_D = 30*np.pi/180
C_phi = 20*np.pi/180

C_lam12 = 1.5
C_lam32 = 0.3
#C_MV_B = 30*np.pi/180

#var_names = ['Time','Half Interval','Bx','By','Bz','Bt','x','y','z','range','tm']
var_names = ['Time','Half Interval','Bx','By','Bz']

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


if csv_file_name[:3]!='DSL':
    print("Ensure B-field data is in DSL coordinates")




################################################################################
    
######################## PROCESSING ############################################

################################################################################


t_overall = df_arr[:,0]
B_x_overall = df_arr[:,2]
B_y_overall = df_arr[:,3]
B_z_overall = df_arr[:,4]
#B_mag = df_arr[:,5]
print("Progress: Variables extracted from dataframe array")

#!!!!!!!!!!!!! artificial offset set here:
if ARTIFICIAL_OFFSET:
    B_z_overall = np.array(len(B_z_overall)*[artificial_Bz_offset]) + np.array(B_z_overall)

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
O_z_overall = []
O_z_overall_overlapping = []

for sheathpass_num in range(0,len(sheath_df_arr)):
    print("Progress:",sheathpass_num)
    
    ts_str = sheath_df_arr[sheathpass_num][0]
    tf_str = sheath_df_arr[sheathpass_num][1]
    
    ts_num = matplotlib.dates.date2num(datetime.strptime(ts_str,'%Y-%m-%dT%H:%M:%S.%fZ'))
    tf_num = matplotlib.dates.date2num(datetime.strptime(tf_str,'%Y-%m-%dT%H:%M:%S.%fZ'))
    sheath_itv_midpoints.append((ts_num+tf_num)/2)


    modified_data = mmcal_functions.data_in_interval([t_overall,B_x_overall,B_y_overall,B_z_overall],[ts_num,tf_num])
    t = modified_data[0][0]
    t_days = modified_data[0][1]
    t_secs = modified_data[0][2]
    B_x = modified_data[1][0]
    B_y = modified_data[1][1]
    B_z = modified_data[1][2]
    B_mag = modified_data[2]
    B_xy =np.array( [(np.sqrt((B_x[i])**2 + (B_y[i])**2 ) ) for i in range(0,len(t))] )
    
        
    
    
    B_xy = np.empty(len(t))
    for i in range(0,len(t)):
        B_xy[i] = np.sqrt((B_x[i])**2 + (B_y[i])**2 )
    
    B_vec = []
    for i in range(0,len(t)):
        B_vec.append(np.array([B_x[i],B_y[i],B_z[i]])/B_mag[i])    
    
    B_polar = mmcal_functions.xyz2polar(B_x,B_y,B_z)
    
    theta_B = B_polar[1]
    phi_B = B_polar[2]
    
    
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
    
    
    
    


    Mmode_count_PN16 = 0
    Mmode_count_SLD08 = 0
    
    Mmode_indices_PN16 = []
    Mmode_indices_SLD08 = []
    Mmode_indices_overlapping = []
    Mmode_indices_PN16_only = []
    Mmode_indices_SLD08_only = []
    
    O_z_PN16 = []

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
            O_z_PN16.append(O_z_unfiltered[i])
            if SLD08_satisfy:
                if len(Mmode_indices_PN16)>0 and len(Mmode_indices_SLD08)>0:
                    if Mmode_indices_PN16[-1]==Mmode_indices_SLD08[-1]:
                        Mmode_indices_overlapping.append(i)           
            else:
                Mmode_indices_PN16_only.append(i)
        else:
            if SLD08_satisfy:
                Mmode_indices_SLD08_only.append(i)
    
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






    count_stats[0].append(Mmode_count_PN16)
    count_stats[1].append(Mmode_count_SLD08)
    count_stats[2].append(len(Mmode_indices_overlapping))
    count_stats[3].append(len(Mmode_sitv_times_PN16_only))
    count_stats[4].append(len(Mmode_sitv_times_SLD08_only))
        

    
    #if len(O_z) > 0:
    if len(Mmode_indices_overlapping) > 0:
        
        sheath_itv_midpoints_plot1.append((ts_num+tf_num)/2)
        
        #Obtaining the pdf using kde
        O_z_PN16_T = np.transpose(np.array([O_z_PN16]))
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
        
        

################################################################################
    
######################## PLOTTING ##############################################

################################################################################

count_stats=np.array(count_stats)
kde_stats=np.array(kde_stats)

kde_overall= KernelDensity(kernel='gaussian',bandwidth=1.0).fit(np.transpose(np.array([O_z_overall])))
log_pdf_overall = kde_overall.score_samples(np.transpose(np.array([O_z_overall])))
pdf_overall = np.exp(log_pdf_overall)

if PLOT:    
    f1=plt.figure()
    ax11 = f1.add_subplot(211)
    ax12 = f1.add_subplot(212)    
    ax11.plot_date(sheath_itv_midpoints_plot1,kde_stats[0][0],markersize=2.0)
    ax12.plot_date(sheath_itv_midpoints_plot1,kde_stats[0][1],markersize=2.0)
    #ax11.set_title(r"$O_{zf}$ estimates per sheath-pass C3 "+plot_time_text+", PN16, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(t_int,shift,artificial_Bz_offset))
    ax11.set_title(r"$O_{zf}$ estimates per day C3 "+plot_time_text+", PN16, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(t_int,shift,artificial_Bz_offset))    
    ax11.set_ylabel("KDE dist maximum")
    ax12.set_ylabel("KDE dist mean")
    
    f2=plt.figure()
    ax21 = f2.add_subplot(211)
    ax22 = f2.add_subplot(212)
    ax21.plot_date(sheath_itv_midpoints_plot1,kde_stats[1][0],markersize=2.0)
    ax22.plot_date(sheath_itv_midpoints_plot1,kde_stats[1][1],markersize=2.0)    
    #ax21.set_title(r"$O_{zf}$ estimates per sheath-pass C3 "+plot_time_text+", PN16&2eigen-conditions, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(t_int,shift,artificial_Bz_offset))
    ax21.set_title(r"$O_{zf}$ estimates per day C3 "+plot_time_text+", PN16&2eigen-conditions, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(t_int,shift,artificial_Bz_offset))    
    ax21.set_ylabel("KDE dist maximum")
    ax22.set_ylabel("KDE dist mean")


    f3=plt.figure()
    ax3 = f3.add_subplot(111)
    ax3.scatter(O_z_overall,pdf_overall,s=1.0,color='orange',label="PN16")
    ax3.set_title("PDF from KDE: C3 "+plot_time_text+", %d PN16 intervals, $(t_{si}=%d,t_{sh}=%d)$ %.2fnT AO"%(len(O_z_overall),t_int,shift,artificial_Bz_offset))
    ax3.set_ylabel("Prob. Density")
    ax3.set_xlabel(r"$O_{z}$ (nT)")
    ax3.legend()        
    
    plt.show()        






