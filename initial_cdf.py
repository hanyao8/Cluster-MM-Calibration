# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 15:03:24 2018

@author: Choon
"""

#some initial cdf reading

import cdflib
import os
import matplotlib.pyplot as plt
import numpy as np

C_xy = 0.3
C_B = 30*np.pi/180
C_D = 30*np.pi/180
C_phi = 20*np.pi/180

#cdf_file1 = cdflib.CDF(os.getcwd()+"\\tha_l1_fgm_20180106_v01.cdf") #06/01/2018 FGM level 1
#cdf_file2 = cdflib.CDF(os.getcwd()+"\\tha_l2_fgm_20180106_v01.cdf") #06/01/2018 FGM level 2

cdf_file1 = cdflib.CDF(os.getcwd()+"\\thc_l1_fgm_20080702_v01.cdf") #02/07/2008 FGM level 1
cdf_file2 = cdflib.CDF(os.getcwd()+"\\thc_l2_fgm_20080702_v01.cdf") #02/07/2008 FGM level 2

#print(cdf_file1.cdf_info())
#print(cdf_file2.cdf_info())
#print(cdf_file2.varinq('tha_fgm_fgs_quality'))
#print(cdf_file2.attinq())

#x1=cdf_file2.varget('tha_fgm_fgs_quality')
fgs_Btot = cdf_file2.varget('thc_fgs_btotal')
fgl_Btot = cdf_file2.varget('thc_fgl_btotal')
fgh_Btot = cdf_file2.varget('thc_fgh_btotal')
fgs_dsl = cdf_file2.varget('thc_fgs_dsl')

fgs_t=cdf_file2.varget('thc_fgs_time')
fgl_t=cdf_file2.varget('thc_fgl_time')
fgh_t=cdf_file2.varget('thc_fgh_time')

fgs_t_itv = fgs_t[int(len(fgs_t)/24*14):int(len(fgs_t)/24*19)]
fgs_dsl_itv = fgs_dsl[int(len(fgs_t)/24*14):int(len(fgs_t)/24*19)]
fgs_Btot_itv = fgs_Btot[int(len(fgs_t)/24*14):int(len(fgs_t)/24*19)]

fgs_Bxy_itv = np.empty(len(fgs_t_itv))
for i in range(0,len(fgs_t_itv)):
    fgs_Bxy_itv[i] = np.sqrt((fgs_dsl_itv[i][0])**2 + (fgs_dsl_itv[i][1])**2 )

t0_itv = fgs_t_itv[0]

subintervals = np.empty(( int((4*3600-170)/10) ,2))
for i in range(0,int((4*3600-170)/10)):
    subintervals[i][0] = t0_itv + i*10
    subintervals[i][1] = t0_itv + i*10 + 180


fgs_Bxy_sitv = []
fgs_Bx_sitv = []
fgs_By_sitv = []
fgs_Bz_sitv = []

fgs_Bx_sitv_mean = np.empty(int((4*3600-170)/10))
fgs_By_sitv_mean = np.empty(int((4*3600-170)/10))
fgs_Bz_sitv_mean = np.empty(int((4*3600-170)/10))

fgs_Bxy_sitv_min = np.empty(int((4*3600-170)/10))
fgs_Bxy_sitv_max = np.empty(int((4*3600-170)/10))
fgs_Bxy_sitv_mean = np.empty(int((4*3600-170)/10))

for i in range(0,int((4*3600-170)/10)):
    fgs_Bx_sitv.append(fgs_dsl_itv[np.argmax(fgs_t_itv>subintervals[i][0]):np.argmax(fgs_t_itv>subintervals[i][1]),0])
    fgs_By_sitv.append(fgs_dsl_itv[np.argmax(fgs_t_itv>subintervals[i][0]):np.argmax(fgs_t_itv>subintervals[i][1]),1])
    fgs_Bz_sitv.append(fgs_dsl_itv[np.argmax(fgs_t_itv>subintervals[i][0]):np.argmax(fgs_t_itv>subintervals[i][1]),2])
    
    fgs_Bxy_sitv.append(fgs_Bxy_itv[np.argmax(fgs_t_itv>subintervals[i][0]):np.argmax(fgs_t_itv>subintervals[i][1])])
    
    
    fgs_Bx_sitv_mean[i] = np.mean(fgs_Bx_sitv[i])
    fgs_By_sitv_mean[i] = np.mean(fgs_By_sitv[i])
    fgs_Bz_sitv_mean[i] = np.mean(fgs_Bz_sitv[i])    
    
    
    fgs_Bxy_sitv_min[i] = min(fgs_Bxy_sitv[i])
    fgs_Bxy_sitv_max[i] = max(fgs_Bxy_sitv[i])
    fgs_Bxy_sitv_mean[i] = np.mean(fgs_Bxy_sitv[i])

#fgs_Bxy_sitv_mean = np.sqrt(fgs_Bx_sitv_mean**2 + fgs_By_sitv_mean**2)






Bxy_fluct = (fgs_Bxy_sitv_max-fgs_Bxy_sitv_min)/fgs_Bxy_sitv_mean


Mmode_count = 0

for i in range(0,len(Bxy_fluct)):
    if Bxy_fluct[i]>C_xy:
        Mmode_count+=1
#to insert some more conditionals here...
#i.e. those associated with phi, theta_B and theta_D.
#To do: Principal Component Analysis


"""
print(len(x2))
print(len(x3))
print(len(x4))
"""
"""
print(x2)
print(x3)
print(x4)
"""

f1=plt.figure()
f2=plt.figure()
f3=plt.figure()

ax1=f1.add_subplot(111)
ax2=f2.add_subplot(111)
ax3=f3.add_subplot(111)

ax1.scatter(fgs_t,fgs_Btot,s=1.0)
ax2.scatter(fgl_t,fgl_Btot,s=1.0)
ax3.scatter(fgh_t,fgh_Btot,s=1.0)
plt.show()

