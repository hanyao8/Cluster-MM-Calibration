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

for jul_day in range(1,32):
    
    if jul_day<10:
        jul_day_str = '0'+str(jul_day)
    else:
        jul_day_str = str(jul_day)

    cdf_file_name = 'thc_l2_fgm_200807'+jul_day_str+'_v01.cdf'
    cdf_file1 = cdflib.CDF(os.getcwd()+'\\cdf_data\\'+cdf_file_name) #02/07/2008 FGM level 1
    
    
    
    
    
    ####################################################################
    
    
    
    
    epoch1970 = 719163.0
    epoch20080702 = 733225.0
    
    fgs_Btot = cdf_file1.varget('thc_fgs_btotal')
    fgs_dsl = cdf_file1.varget('thc_fgs_dsl')
    fgs_t=cdf_file1.varget('thc_fgs_time')
    #fgs_t2 = (fgs_t - np.array([fgs_t[0]]*len(fgs_t))) /3600/24
    fgs_t2 = fgs_t/3600/24 + np.array([epoch1970]*len(fgs_t) ) 
    
    
    #matplotlib.dates.num2date(1214957297.0989337/3600/24 + epoch1970)
    
    #t = fgs_t2 + np.array([epoch20080702]*len(fgs_t))
    t_pltobj = a=matplotlib.dates.num2date(fgs_t2)
    t_str = [x.strftime('%Y-%m-%dT%H:%M:%S.%fZ') for x in t_pltobj]
    B_x = fgs_dsl[:,0]
    B_y = fgs_dsl[:,1]
    B_z = fgs_dsl[:,2]
    B_mag = fgs_Btot
    
    
    file_name = "DSL_THEMIS_C_FGM_SPINRES_2008Jul"
    csv_file = io.open(os.getcwd() +"\\"+ file_name + ".csv",'a')
    
    for i in range(0,len(t_str)):
        data_row_str = t_str[i] +",0,"+ str(B_x[i]) +","+ str(B_y[i]) +","+ str(B_z[i])
        csv_file.write(data_row_str+"\n")        
    
    




