import numpy as np
import matplotlib
from datetime import datetime

#import pywt
#from scipy.integrate import quad

#from sklearn.neighbors.kde import KernelDensity


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
    
def data_in_interval(data,se_times):
    t = data[0]
    B_x = data[1]    
    B_y = data[2]
    B_z = data[3]
    
    data_start_time = se_times[0]
    data_end_time = se_times[1]
    
    t_datetime = []
    for i in range(0,len(t)):
        strpdtime=datetime.strptime(t[i],'%Y-%m-%dT%H:%M:%S.%fZ')
        t_datetime.append(strpdtime)
    t_days = matplotlib.dates.date2num(t_datetime)
    t_secs= t_days*24*3600


    if data_start_time > t_days[0]:
        data_start_index = np.argmax(t_days>data_start_time)
    else:
        data_start_index = 0
    if data_end_time < t_days[-1]:
        data_end_index = np.argmax(t_days>data_end_time)
    else:
        data_end_index = len(t_days)             
    
    t = t[data_start_index:data_end_index]
    B_x = B_x[data_start_index:data_end_index]
    B_y = B_y[data_start_index:data_end_index]
    B_z = B_z[data_start_index:data_end_index]
    #B_mag = B_mag[data_start_index:data_end_index]
    B_x.astype(float)
    B_y.astype(float)
    B_z.astype(float)
    B_mag = (B_x**2+B_y**2+B_z**2)**(0.5)
    
    t_days = t_days[data_start_index:data_end_index]
    t_secs = t_secs[data_start_index:data_end_index]
    
    return([[t,t_days,t_secs],[B_x,B_y,B_z],B_mag])
    
    
    
def sitvfy(input_data,settings_params):
    t_secs = input_data[0]
    B_x = input_data[1][0]
    B_y = input_data[1][1]
    B_z = input_data[1][2]
    B_mag = input_data[2]    
    
    t_int = settings_params[0]
    shift = settings_params[1]
    
    B_xy = []
    for i in range(0,len(t_secs)):
        B_xy.append (np.sqrt((B_x[i])**2 + (B_y[i])**2 ) )
    B_xy = np.array(B_xy)
    
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
    Bmag_sitv_mean = []
    
    Bxy_sitv_min = []
    Bxy_sitv_max = []
    Bxy_sitv_mean = []
    
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
            Bmag_sitv_mean.append(np.mean(B_mag[si_start:si_end]))
    
            Bxy_sitv_min.append(min(Bxy_sitv[-1]))
            Bxy_sitv_max.append(max(Bxy_sitv[-1]))
            Bxy_sitv_mean.append(np.mean(Bxy_sitv[-1]))
    
    Bx_sitv_mean = np.array(Bx_sitv_mean)
    By_sitv_mean = np.array(By_sitv_mean)
    Bz_sitv_mean = np.array(Bz_sitv_mean)             
            
    subintervals = np.array(gapadj_subintervals.copy())
    sitv_midpoints_secs = (subintervals[:,0]+subintervals[:,1])/2
    
    return([sitv_midpoints_secs,[Bx_sitv,By_sitv,Bz_sitv],Bmag_sitv_mean,[Bxy_sitv_min,Bxy_sitv_max,Bxy_sitv_mean]])
    
    

def MVA(input_data):
    
    sitv_midpoints_secs = input_data[0]
    Bx_sitv = input_data[1][0]
    By_sitv = input_data[1][1]    
    Bz_sitv = input_data[1][2]
    Bxy_sitv_min = input_data[2][0]
    Bxy_sitv_max = input_data[2][1]   
    Bxy_sitv_mean = input_data[2][2]    
    Bx_sitv_mean = [np.mean(Bx_sitv[i]) for i in range(0,len(sitv_midpoints_secs))]
    By_sitv_mean = [np.mean(By_sitv[i]) for i in range(0,len(sitv_midpoints_secs))]
    Bz_sitv_mean = [np.mean(Bz_sitv[i]) for i in range(0,len(sitv_midpoints_secs))]

    Bxy_fluct_PN16 = (np.array(Bxy_sitv_max)-np.array(Bxy_sitv_min))/np.array(Bxy_sitv_mean)
    phi_PN16 = []
    theta_B_PN16 = []
    theta_D_PN16 = []
    
    B_x1_angle = []
    theta_D = []
    phi_D = []
    lam1_lam2 = []
    lam3_lam2 = []
    
    O_z_unfiltered = []
    
    for i in range(0,len(sitv_midpoints_secs)):
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
        
        lam1_lam2.append(lam1/lam2)
        lam3_lam2.append(lam3/lam2)
    
        x1 = eigen[1][:,lam1_index]
        x1_xy = np.sqrt(x1[0]**2+x1[1]**2)
    
        B_dir = np.array([Bx_sitv_mean[i],By_sitv_mean[i],Bz_sitv_mean[i]])
        B_dir /= np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2)
        B_dir_xy = np.sqrt(B_dir[0]**2+B_dir[1]**2)
    
        if np.dot(x1,B_dir) < 0:
            x1=-x1
    
        theta_D.append( (np.arctan2( [x1_xy],[x1[2]] ))[0])
        phi_D.append( (np.arctan2( [x1[1]],[x1[0] ]) )[0] )
        
    
        B_x1_angle.append( np.arccos(np.dot(x1,B_dir)))
        phi_PN16.append( np.arccos(np.dot([x1[0],x1[1]],[B_dir[0],B_dir[1]]) \
                                         /np.sqrt(x1[0]**2+x1[1]**2)/np.sqrt(B_dir[0]**2+B_dir[1]**2) ))
        """           
        else:
            B_x1_angle.append( np.arccos(np.dot(x1,B_dir))  )
            phi_PN16.append( np.arccos(np.dot([x1[0],x1[1]],[B_dir[0],B_dir[1]]) \
                                             /np.sqrt(x1[0]**2+x1[1]**2)/np.sqrt(B_dir[0]**2+B_dir[1]**2) ) )
        """
    
        theta_B_PN16.append( np.arctan2(B_dir[2],B_dir_xy) )
        theta_D_PN16.append( np.arctan2(x1[2],x1_xy) )
        
        O_z_unfiltered.append( Bz_sitv_mean[i] - x1[2]/x1_xy*Bxy_sitv_mean[i] )
                        
    return([[phi_PN16,theta_D_PN16,theta_B_PN16,Bxy_fluct_PN16]\
                ,[B_x1_angle,theta_D,phi_D,lam1_lam2,lam3_lam2]\
                ,O_z_unfiltered] )
    
    
    
    
def MVA_offsetref(input_data):
    
    sitv_midpoints_secs = input_data[0]
    Bx_sitv = input_data[1][0]
    By_sitv = input_data[1][1]    
    Bz_sitv = input_data[1][2]
    Bxy_sitv_min = input_data[2][0]
    Bxy_sitv_max = input_data[2][1]   
    Bxy_sitv_mean = input_data[2][2]    
    Bx_sitv_mean = [np.mean(Bx_sitv[i]) for i in range(0,len(sitv_midpoints_secs))]
    By_sitv_mean = [np.mean(By_sitv[i]) for i in range(0,len(sitv_midpoints_secs))]
    Bz_sitv_mean = [np.mean(Bz_sitv[i]) for i in range(0,len(sitv_midpoints_secs))]

    Bxy_fluct_PN16 = (np.array(Bxy_sitv_max)-np.array(Bxy_sitv_min))/np.array(Bxy_sitv_mean)
    phi_PN16 = []
    theta_B_PN16 = []
    theta_D_PN16 = []
    
    B_x1_angle = []
    theta_D = []
    phi_D = []
    lam1_lam2 = []
    lam3_lam2 = []

    
    O_z_unfiltered = []
    dO=[]
    
    d_g=10**(-4)
    d_bn=10**(-12)
    d_thetaD=[]    
    d_b=[]
    bxymirror=[]
    bzmirror=[]
    d_thetaB=[]
    theta_mb=[]
    theta_md=[]
    test=[]
    oz1=[]    
    
    for i in range(0,len(sitv_midpoints_secs)):
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
        
        lam1_lam2.append(lam1/lam2)
        lam3_lam2.append(lam3/lam2)
    
        x1 = eigen[1][:,lam1_index]
        x1_xy = np.sqrt(x1[0]**2+x1[1]**2)
    
        B_dir = np.array([Bx_sitv_mean[i],By_sitv_mean[i],Bz_sitv_mean[i]])
        B_dir /= np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2)
        B_dir_xy = np.sqrt(B_dir[0]**2+B_dir[1]**2)
    
        if np.dot(x1,B_dir) < 0:
            x1=-x1
    
        theta_D.append( (np.arctan2( [x1_xy],[x1[2]] ))[0])
        phi_D.append( (np.arctan2( [x1[1]],[x1[0] ]) )[0] )
        B_magn=np.sqrt(Bx_sitv_mean[i]**2+By_sitv_mean[i]**2+Bz_sitv_mean[i]**2)        
    
        B_x1_angle.append( np.arccos(np.dot(x1,B_dir)))
        phi_PN16.append( np.arccos(np.dot([x1[0],x1[1]],[B_dir[0],B_dir[1]]) \
                                         /np.sqrt(x1[0]**2+x1[1]**2)/np.sqrt(B_dir[0]**2+B_dir[1]**2) ))
        """           
        else:
            B_x1_angle.append( np.arccos(np.dot(x1,B_dir))  )
            phi_PN16.append( np.arccos(np.dot([x1[0],x1[1]],[B_dir[0],B_dir[1]]) \
                                             /np.sqrt(x1[0]**2+x1[1]**2)/np.sqrt(B_dir[0]**2+B_dir[1]**2) ) )
        """
    
        theta_B_PN16.append( np.arctan2(B_dir[2],B_dir_xy) )
        theta_D_PN16.append( np.arctan2(x1[2],x1_xy) )
        
        O_z_unfiltered.append( Bz_sitv_mean[i] - x1[2]/x1_xy*Bxy_sitv_mean[i] )
                        
    return([[phi_PN16,theta_D_PN16,theta_B_PN16,Bxy_fluct_PN16]\
                ,[B_x1_angle,theta_D,phi_D,lam1_lam2,lam3_lam2]\
                ,O_z_unfiltered] )
    

def get_Mmode_times(sitv_midpoints_days,Mmode_indices_PN16_only\
                                     ,Mmode_indices_SLD08_only\
                                     ,Mmode_indices_overlapping):
    
    #Mmode_sitv_times_PN16 = []
    #Mmode_sitv_times_SLD08 = []
    Mmode_sitv_times_PN16_only = []
    #Mmode_sitv_times_PN16_days = []
    Mmode_sitv_times_SLD08_only = []
    #Mmode_sitv_times_SLD08_days = []    
    Mmode_sitv_times_overlapping = [] 
    #Mmode_sitv_times_overlapping_days = [] 
    
    #for i in range(0,len(Mmode_indices_PN16)):
        #Mmode_sitv_times_PN16.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_PN16[i]]))
        #Mmode_sitv_times_PN16_days.append(sitv_midpoints_days[Mmode_indices_PN16[i]])
    #for i in range(0,len(Mmode_indices_SLD08)):
        #Mmode_sitv_times_SLD08.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_SLD08[i]]))   
        #Mmode_sitv_times_SLD08_days.append(sitv_midpoints_days[Mmode_indices_SLD08[i]])
    for i in range(0,len(Mmode_indices_PN16_only)):
        Mmode_sitv_times_PN16_only.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_PN16_only[i]]))
    for i in range(0,len(Mmode_indices_SLD08_only)):
        Mmode_sitv_times_SLD08_only.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_SLD08_only[i]]))   
    for i in range(0,len(Mmode_indices_overlapping)):        
        Mmode_sitv_times_overlapping.append(matplotlib.dates.num2date(sitv_midpoints_days[Mmode_indices_overlapping[i]]))
        #Mmode_sitv_times_overlapping_days.append(sitv_midpoints_days[Mmode_indices_overlapping[i]])
    
    return([Mmode_sitv_times_PN16_only,Mmode_sitv_times_SLD08_only,Mmode_sitv_times_overlapping])
    


def sitv_in_sheath_frac(sitv_midpoints_days,indices_arrs,sheath_times):

    Mmode_indices_PN16 = indices_arrs[0]
    Mmode_indices_SLD08 = indices_arrs[1]
    Mmode_indices_overlapping = indices_arrs[2]
    
    sheath_start_time = sheath_times[0]
    sheath_end_time = sheath_times[1]
    
    Mmode_sitv_times_PN16_days = [(sitv_midpoints_days[Mmode_indices_PN16[i]]) for i in range(0,len(Mmode_indices_PN16))]
    Mmode_sitv_times_SLD08_days = [(sitv_midpoints_days[Mmode_indices_SLD08[i]]) for i in range(0,len(Mmode_indices_SLD08))]  
    Mmode_sitv_times_overlapping_days =[(sitv_midpoints_days[Mmode_indices_overlapping[i]]) for i in range(0,len(Mmode_indices_overlapping))]
    
    PN16_frac_in_sheath = (np.argmax(np.array(Mmode_sitv_times_PN16_days)>sheath_end_time) - np.argmax(np.array(Mmode_sitv_times_PN16_days)>sheath_start_time) )/len(Mmode_indices_PN16)
    SLD08_frac_in_sheath = (np.argmax(np.array(Mmode_sitv_times_SLD08_days)>sheath_end_time) - np.argmax(np.array(Mmode_sitv_times_SLD08_days)>sheath_start_time) )/len(Mmode_indices_SLD08)
    overlapping_frac_in_sheath = (np.argmax(np.array(Mmode_sitv_times_overlapping_days)>sheath_end_time) - np.argmax(np.array(Mmode_sitv_times_overlapping_days)>sheath_start_time) )/len(Mmode_indices_overlapping)
    
    print(PN16_frac_in_sheath)
    print(SLD08_frac_in_sheath)
    print(overlapping_frac_in_sheath)
    
    return()
	
	
def dailyOff(O_z,starttim):
    
    #takes in a list of O_z and separate them into days, returns in list form. Together
    #with the associated error
    O_z = np.array(O_z)
    O_zd=[]
    errord=[]
    for i in range(len(np.unique(starttim))):
        O_zd.append([])
        errord.append([])
    
    count=0
    for i in range(len(O_z)):
            if i==0 and starttim[0]==starttim[1]:
                O_zd[count].append(O_z[i])
                errord[count].append(error[i])
            elif starttim[i-1]==starttim[i]:
                O_zd[count].append(O_z[i])
                errord[count].append(error[i])
            elif count<len(O_zd)-1:
                count+=1
                
                
    for i in errord:
        i=np.array(i)  
    for i in O_zd:
        i=np.array(i)    
    
    O_ztd=[]
    
    for i in range(len(O_zd)):
        O_ztd.append(np.transpose(np.array(np.array(O_zd[i]))))
    
    return O_ztd,errord
	
def kde(mean,sigma,weight):    
    x=np.linspace(-15,15,1000)
    gaussian=[]
    ker=0
    for i in range(len(mean)):
        mu, width=mean, sigma
        y=1/(width[i] * np.sqrt(2 * np.pi)*weight[i])*np.exp( - (x - mu[i])**2 / (2 * width[i]**2))
        gaussian.append(y)
        ker+=y
    ker=ker/len(mean)
    plt.plot(x,ker,color='blue',label='weighted kde h=1')
    plt.legend()
    return x[np.where(ker==np.max(ker))[0][0]]
	
def adapkde(mean):
    h=bandwidth.sskernel(mean)[2]
    x=np.linspace(-20,20,100000)
    ker=0
    for i in range(len(mean)):
        mu, width=mean, sigma
        y=1/(h * np.sqrt(2 * np.pi))*np.exp( - (x - mu[i])**2 / (2 * h**2))
        ker+=y
    ker=ker/len(mean)
    plt.plot(x,ker,color='green',label='adaptive kde no weighting')
    plt.legend()
    plt.title("Types of KDE")
    plt.xlabel("$O_z/nT$")
    plt.ylabel("P/$nT^{-1}$")
    plt.show()    
    return x[np.where(ker==np.max(ker))[0][0]]
    
def wadapkde(mean,weight):
    h=bandwidth.sskernel(mean)[2]
    x=np.linspace(-20,20,100000)
    ker=0
    for i in range(len(mean)):
        mu, width=mean, sigma
        y=1/(h * np.sqrt(2 * np.pi)*weight[i])*np.exp( - (x - mu[i])**2 / (2 * h**2))
        ker+=y
    ker=ker/len(mean)
    """
    plt.plot(x,ker,color='red',label='adaptive kde with weighting')
    plt.legend()
    plt.show()  
    """
    return x[np.where(ker==np.max(ker))[0][0]]


	
