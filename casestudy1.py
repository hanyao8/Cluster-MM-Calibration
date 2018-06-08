#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 14:19:09 2018

@author: jia_qu
"""

gamma=np.pi/5
k=(2*np.pi)/10
lamb=10
w_0=(1/10.0)*lamb
X=np.linspace(-4000,4000,100)
Y=np.linspace(-4000,4000,100)
u_0=1
rho_0=400*lamb

normIx=[]
normExxyy=[]
normExxxx=[]




for i in np.linspace(0,100*lamb,10000):
    
        xi=np.sqrt(1+((4*(i**2))/((k**2)*(w_0**4))))
        
        w=w_0*xi
        C=(u_0*(w_0/w)*np.cos(gamma))/np.sqrt(2)
        r=(np.sqrt(2)*np.tan(gamma)*rho)/w
        
        phi_x=np.arctan2(2*i,(k*w_0**2))
        
        norm=1+r**2
        renorm=(np.pi*u_0**2*w_0**2*np.cos(gamma)**2)*0.25*(1+np.tan(gamma)**2)
        
        I_5=(w/2)**6*(1-0.5*(((2*rho_0/w)**4)+2*((2*rho_0/w)**3)+2)*np.exp((-4*rho_0**2)/w**2))
        
     
        
        I_xxxx=(w**4*(13+4*np.cos(2*gamma)-np.cos(4*gamma))+np.exp((-4*rho_0**2)/w**2)*(w**4*
                (-13-4*np.cos(2*gamma)+np.cos(4*gamma))-16*rho_0**2*w**2*(5+3*np.cos(2*gamma))*np.sin(gamma)**2
                -64*rho_0**4*(np.sin(gamma))**4))
        

        
        
        
        #g_xx
        Itrial=np.cos(gamma)**2*(1-np.exp(-2*rho_0**2)/(w**2)+(np.tan
            (gamma)**2/w**2)*(w**2-w**2*np.exp((-2*rho_0**2)/(w**2))-2*rho_0**2*np.exp((-2*rho_0**2)/(w**2))))
        normIx.append(Itrial)
        
        #g_xxyy
        E_xxyy=(8*np.cos(gamma)**4)/(np.pi*w**4)*((w**2/4)*(1-np.exp((-2*rho_0**2)/(w**2)))+(4/w**4)*np.tan(gamma)**4*I_5)
        normExxyy.append(E_xxyy)
        
        #g_xxxx
        f=np.exp((4*rho_0**2)/(w**2))-1
        E_xxxx=(np.cos(gamma)**4/(2*(w**6)*np.pi))*np.exp((-4*rho_0**2)/(w**2))*(
                2*(f)*w**4-4*w**2*(4*rho_0**2-(f)*w**2)*np.tan(gamma)**2-(8*rho_0**4
                  -4*rho_0**2*w**2-f*w**4)*np.tan(gamma)**4)
        
        normExxxx.append(E_xxxx)