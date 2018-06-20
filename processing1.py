# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 15:57:17 2018

@author: Choon
"""

#processing

import numpy as np
import matplotlib
from datetime import datetime
import pywt
from scipy.integrate import quad

########################## FUNCTIONS #####################################

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

########################## FUNCTIONS #####################################



class processing:
    def __init__(self,df_arr,params,SETTING):

        t = df_arr[:,0]
        B_x = df_arr[:,2]
        B_y = df_arr[:,3]
        B_z = df_arr[:,4]
        B_mag = df_arr[:,5]

        data_start_time = params[0][0]
        data_end_time = params[0][1]
        
        t_int = params[1][0]
        shift = params[1][1]
        
        if SETTING=='themis_validation':
            t_days = t
        
        else:
            t_datetime = []
            for i in range(0,len(t)):
                strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
                t_datetime.append(strpdtime)
            t_days = matplotlib.dates.date2num(t_datetime)
        
        t_secs= t_days*24*3600

        data_start_index = np.argmax(t_days>data_start_time)
        data_end_index = np.argmax(t_days>data_end_time)-1
             
        
        t = t[data_start_index:data_end_index]
        B_x = B_x[data_start_index:data_end_index]
        B_y = B_y[data_start_index:data_end_index]
        B_z = B_z[data_start_index:data_end_index]
        B_mag = B_mag[data_start_index:data_end_index]
        

        print(t[0])
        self.t=t

        self.B_x=B_x
        self.B_y=B_y
        self.B_z=B_z
        self.B_mag=B_mag

        #t_datetime = []
        #for i in range(0,len(t)):
        #    strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
        #    t_datetime.append(strpdtime)
        #t_days = matplotlib.dates.date2num(t_datetime)
        #t_secs= t_days*24*3600
        
        t_days = t_days[data_start_index:data_end_index]
        t_secs = t_secs[data_start_index:data_end_index]
        
        self.t_days=t_days
        self.t_secs=t_secs
        
        #print(t_days[0])
        #print(self.t_days[0])

        
        B_xy = np.empty(len(t))
        for i in range(0,len(t)):
            B_xy[i] = np.sqrt((B_x[i])**2 + (B_y[i])**2 )

        B_vec = []
        for i in range(0,len(t)):
            B_vec.append(np.array([B_x[i],B_y[i],B_z[i]])/B_mag[i])    
    
        B_polar = xyz2polar(self.B_x,self.B_y,self.B_z)
        
        self.theta_B = B_polar[1]
        self.phi_B = B_polar[2]




        
        subintervals = np.empty(( int((t_secs[-1]-t_secs[0]-t_int+shift)/shift) ,2))
        for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
            subintervals[i][0] = t_secs[0] + i*shift
            subintervals[i][1] = t_secs[0] + i*shift + t_int

        
        Bx_sitv = []
        By_sitv = []
        Bz_sitv = []
        Bmag_sitv=[]
        Bxy_sitv = []

        self.Bx_sitv_mean = []
        self.By_sitv_mean = []
        self.Bz_sitv_mean = []
        
        Bxy_sitv_min = []
        Bxy_sitv_max = []
        self.Bxy_sitv_mean = []

        gapadj_subintervals = []
        for i in range(0,int((t_secs[-1]-t_secs[0]-t_int+shift)/shift)):
            si_start = np.argmax(t_secs>subintervals[i][0])
            si_end = np.argmax(t_secs>subintervals[i][1])
        
            if len(B_x[si_start:si_end])!=0:
                gapadj_subintervals.append(subintervals[i])
                
                Bx_sitv.append(self.B_x[si_start:si_end])
                By_sitv.append(self.B_y[si_start:si_end])
                Bz_sitv.append(self.B_z[si_start:si_end])
                Bmag_sitv.append(self.B_mag[si_start:si_end])
                Bxy_sitv.append(B_xy[si_start:si_end])
                
                self.Bx_sitv_mean.append(np.mean(Bx_sitv[-1]))
                self.By_sitv_mean.append(np.mean(By_sitv[-1]))
                self.Bz_sitv_mean.append(np.mean(Bz_sitv[-1]))
                
                Bxy_sitv_min.append(min(Bxy_sitv[-1]))
                Bxy_sitv_max.append(max(Bxy_sitv[-1]))
                self.Bxy_sitv_mean.append(np.mean(Bxy_sitv[-1]))
                
                
        subintervals = np.array(gapadj_subintervals.copy())
        self.subintervals = subintervals


        self.B_x1_angle = []
        self.theta_D = []
        self.phi_D = []
        self.lam1_lam2 = []
        self.lam3_lam2 = []
        
        self.Bxy_fluct_PN16 = (np.array(Bxy_sitv_max)-np.array(Bxy_sitv_min))/self.Bxy_sitv_mean
        self.phi_PN16 = []
        self.theta_B_PN16 = []
        self.theta_D_PN16 = []
        
        for i in range(0,len(subintervals)):
            data = np.array([Bx_sitv[i],By_sitv[i],Bz_sitv[i]])
            eigen = np.linalg.eig(varlist(data))   
    
            lam1_index = np.argmax(eigen[0])
            lam3_index = np.argmin(eigen[0])
            for i in [0,1,2]:
                if i!=lam1_index and i!=lam3_index:
                    lam2_index=i
                    break
                
            lam1 = eigen[0][lam1_index]
            lam2 = eigen[0][lam2_index]
            lam3 = eigen[0][lam3_index]
            
            self.lam1_lam2.append(lam1/lam2)
            self.lam3_lam2.append(lam3/lam2)
    
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

    
            B_dir = np.array([self.Bx_sitv_mean[i],self.By_sitv_mean[i],self.Bz_sitv_mean[i]])
            B_dir /= np.sqrt(self.Bx_sitv_mean[i]**2+self.By_sitv_mean[i]**2+self.Bz_sitv_mean[i]**2)
            B_dir_xy = np.sqrt(B_dir[0]**2+B_dir[1]**2)
            
            self.B_x1_angle.append( np.arccos(np.dot(x1,B_dir)) * 180/np.pi )
            #self.B_x1_angle.append( np.arccos(abs(np.dot(x1,B_dir))) * 180/np.pi )

            self.theta_D.append( (np.arctan2( [x1_xy],[x1[2]] ))[0])
            self.phi_D.append( (np.arctan2( [x1[1]],[x1[0] ]) )[0] )
            
            self.phi_PN16.append( np.arccos(np.dot([x1[0],x1[1]],[B_dir[0],B_dir[1]]) \
                                             /np.sqrt(x1[0]**2+x1[1]**2)/np.sqrt(B_dir[0]**2+B_dir[1]**2) ) )
            
            self.theta_B_PN16.append( np.arctan2(B_dir[2],B_dir_xy) )
            self.theta_D_PN16.append( np.arctan2(x1[2],x1_xy) )          
                
        if SETTING=='peakness':
            widths = np.linspace(1, 61,100)
            
            (cA, cD) = pywt.cwt(B_mag, widths,wavelet='mexh')

            b=[]
            bsitv=[]
            bsitv_mean = []
            peakness= []         
            #bsitv_mean = np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))
            #peakness= np.empty(int((t_secs[-1]-t_secs[0]-t_int+shift)/shift))   
            deltas=widths[1]-widths[0]
            val=0
            for i in range(len(B_mag)):
                for j in range(len(widths)):
                    val+=cA[j][i]*(widths[j])**(-3/2)*deltas
                    
                b.append(val)
                
            for i in range(0,len(subintervals)):
                si_start = np.argmax(t_secs>subintervals[i][0])
                si_end = np.argmax(t_secs>subintervals[i][1])                
    
                bsitv.append(b[si_start:si_end])
                bsitv_mean.append(np.mean(bsitv[i]))
                peakness.append(np.mean((bsitv[i]-bsitv_mean[i])**3)/(np.mean((bsitv[i]-bsitv_mean[i])**2))**1.50)
                
            self.peakness = np.array(peakness)
            self.b_wav = np.array(b)



