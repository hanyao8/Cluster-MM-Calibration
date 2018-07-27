#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 17:59:24 2018

@author: jia_qu
"""
import matplotlib.pyplot as plt

months=["Jan","Feb","March","April"]
month=["Sept","Oct","Nov","Dec","Jan","Feb","March","April"]
madapk=[1.102,0.170,-0.016]
mkde=[0.991,0.152,0.068]
adapk=[0.0386,1.048,0.586,0.126]
kde=[0.225,1.066,0.629,0.736]
s=["2","3","4","5"]



fig,ax = plt.subplots()
for i, txt in enumerate(s):
    ax.annotate(txt, (months[i],adapk[i]))
ax.plot(months,kde,color='b')
plt.plot(month,mkde,color='y')
ax.scatter(months,adapk,color='g',label="Algorithm adaptive kde")
ax.plot(months,adapk,color='g')
ax.plot(months,kde,color='b')
ax.scatter(months,kde,color='b',label="Algorithm kde")
ax.scatter(month,madapk,color='r',label="Manual adaptive kde")
ax.plot(month,madapk,color='r')
ax.scatter(month,mkde,color='y',label="Manual kde")

ax.set_title("Offset without artifitial shift")
ax.set_xlabel("Months")
ax.set_ylabel("$O_z$")
ax.legend()


offadapk=[4.638,6.095,6.095, 5.069]
offkde=[4.820,6.20,6.201 ,4.279]

plt.figure()
plt.title("Offset with $5nT$ artifitial shift")
plt.xlabel("Months")
plt.ylabel("$O_z/nT$")
plt.scatter(months,offadapk,color='g',label="Algorithm adaptive kde")
plt.plot(months,offadapk,color='g')
plt.plot(months,offkde,color='b')
plt.scatter(months,offkde,color='b',label="Algorithm kde")
plt.legend()


#february applied vs calculated offset
plt.figure()
applied=[0,5,-0.95,-0.9938]
wa=[0.886 , 5.691,0.0726,0.03]
kd=[1.009 ,5.811,0.01500,-0.0238]
lowererrora=[2.9584295842958426,3.4600346003460025,2.45,2.46]
uppererrora=[ 2.21,2.57,2.17, 2.175]
aae=[lowererrora,uppererrora]

lowererror=[3.25,4.15,3.15,3.147]
uppererror=[2.74,2.786,2.712,2.711]
er=[lowererror,uppererror]
plt.figure()
plt.plot(applied,wa,color='b',label="adaptive kde")
plt.errorbar(applied,wa,color='b',yerr=aae,fmt='o')
plt.plot(applied,kd,color='r',label="kde")
plt.errorbar(applied,kd,color='r',yerr=er,fmt='o')
plt.legend()
plt.xlabel("Applied offset$/{nT}$")
plt.ylabel("Calculated offset$/{nT}$")
plt.title("February calibration")

#plot of whole months
madapk1=[-1.779,-0.6606,0.0162,0.2233,-0.0338,0.886,0.141,0.0690	]
mkde1=[-1.2338,-0.927,0.167,0.0906,0.191,1.009,0.2298,0.08420]
lowererrora1=[0,0,0,0,-0.001,2.958,2.94,2.008]
uppererrora1=[0,0,0,0,2.4408,2.21,3.19,2.619]
#error for kde
lowererror1=[0,0,0,0,2.546,3.25,3.59,2.6192]
uppererror1=[0,0,0,0,2.7324,2.74,3.487,2.776]
adaperror=[lowererrora1,uppererrora1]
kderror=[lowererror1,uppererror1]
s1=[1682,2512,10445,17483,5922,8285,6877,5659]
fig,ax = plt.subplots()
for i, txt in enumerate(s1):
    ax.annotate(txt, (month[i],madapk1[i]))
ax.plot(month,mkde1,color='b',label="kde")
ax.plot(month,madapk1,color='r',label="adaptive kde")
#ax.scatter(month,mkde1,color='b',label="kde")
#ax.scatter(month,madapk1,color='r',label="adaptive kde")
ax.errorbar(month,madapk1,color='r',yerr=adaperror,fmt='o')
ax.errorbar(month,mkde1,color='b',yerr=kderror,fmt='o')
ax.set_title("Offset without artifitial shift")
ax.set_xlabel("Months")
ax.set_ylabel("$O_z/{nT}$")
ax.legend()





#plots for April


madapk4=[ 0.0690,0.018200,0.0066]
mkde4=[0.0842,0.0234,0.0074]
applied=[0,-0.076, -0.0766-0.0208]
lowererrora4=[2.008,2.022,2.034]
uppererrora4=[2.4408,2.4208,2.4204]
lowererror4=[2.619,2.626,2.632]
uppererror4=[2.776,2.7729,2.773]
adaperror4=[lowererrora4,uppererrora4]
kderror4=[lowererror4,uppererror4]
s4=[5671,5686,5686]
fig,ax = plt.subplots()
for i, txt in enumerate(s4):
    ax.annotate(txt, (applied[i],madapk4[i]))
ax.plot(applied,madapk4,color='r',label="adaptive kde")
ax.plot(applied,mkde4,color='b',label="kde")
#ax.scatter(month,mkde1,color='b',label="kde")
#ax.scatter(month,madapk1,color='r',label="adaptive kde")
ax.errorbar(applied,madapk4,color='r',yerr=adaperror4,fmt='o')
ax.errorbar(applied,mkde4,color='b',yerr=kderror4,fmt='o')
ax.set_title("April calibration")
ax.set_xlabel("Applied offset$/{nT}$")
ax.set_ylabel("Calculated offset$/{nT}$")
ax.legend()

#January applied vs calculated offset

madapk1=[ -0.0338,-0.1074,-0.1138,-0.115]
mkde1=[0.191,0.1238,0.1158,0.115]
applied=[0,-0.0786, -0.0084,-0.001]
lowererrora1=[-0.001,1.760,1.7640,1.7640]
uppererrora1=[2.4408,2.4208,2.4204,2.411]
lowererror1=[2.546,2.552,2.553,2.553]
uppererror1=[2.7324,2.7280,2.731227,2.731227]
adaperror1=[lowererrora1,uppererrora1]
kderror1=[lowererror1,uppererror1]
s1=[5922,5927,5927,5926]
fig,ax = plt.subplots()
for i, txt in enumerate(s1):
    ax.annotate(txt, (applied[i],madapk1[i]))
ax.plot(applied,madapk1,color='r',label="adaptive kde")
ax.plot(applied,mkde1,color='b',label="kde")
#ax.scatter(month,mkde1,color='b',label="kde")
#ax.scatter(month,madapk1,color='r',label="adaptive kde")
ax.errorbar(applied,madapk1,color='r',yerr=adaperror1,fmt='o')
ax.errorbar(applied,mkde1,color='b',yerr=kderror1,fmt='o')
ax.set_title("January calibration")
ax.set_xlabel("Applied offset$/{nT}$")
ax.set_ylabel("Calculated offset$/{nT}$")
ax.legend()
