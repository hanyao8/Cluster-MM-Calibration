# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 14:29:46 2018

@author: jia_qu
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


data=pd.read_csv("/Users/jia_qu/Desktop/2006March.csv",error_bad_lines=False)
data["newtime"]=pd.Series(1, index=data.index)

"""
plt.plot(data.Time,data.Bt)
plt.xlabel("Time")
plt.ylabel("Magnetic field/$nT$")
"""

def varlist(data):
    """Calculate the mean of the magnetic field components and return the variance matrix as a (3,3) array"""
    
    Bxm=data["Bx"].mean()
    Bym=data["By"].mean()
    Bzm=data["Bz"].mean()
    Bxsqm=(data["Bx"]**2).mean()
    Bysqm=(data["By"]**2).mean()
    Bzsqm=(data["Bz"]**2).mean()
    Bxym=(data["Bx"]*data["By"]).mean()
    Bxzm=(data["Bx"]*data["Bz"]).mean()
    Byzm=(data["By"]*data["Bz"]).mean()
    
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
    
    
    
np.linalg.eig(varlist(data))    

    
    
    

#solve eigenvalue equation