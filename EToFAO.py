# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 17:19:58 2018

@author: Omar
Already tested. Stable version
"""
import numpy as np 

 
def EToFAO(ta,u,hr,Rs,doy,lat,deltat,AneHeight=2):
    """
    ETo FAO56
    """
   
    
    
    
    Rs_day =(np.sum(Rs)*3600*deltat)/1e6#(MJ m-2 day-1) daily solar radiation
    
        
    #compute extraterrestrial radiation
    latr=lat*np.pi/180
    dr=1+0.033*np.cos(2*np.pi*doy/365)
    sdec=(23.45*np.cos((2*np.pi*(doy-172)/365)))*np.pi/180#(degree)solar declination angle
    hs=np.arccos(-np.tan(latr)*np.tan(sdec))#half of day duration
    
    rext=37.4*dr*(np.sin(latr)*np.sin(sdec)*hs+np.cos(latr)*np.cos(sdec)*np.sin(hs))
        
    
    tra=Rs_day/rext
    fhs=(tra-0.25)/0.5
    u_m=np.mean(u)
    
    if AneHeight == 2:
        U2=u_m
    else:
        U2=u_m*(4.87/(np.log(67.8*AneHeight-5.42)))
        
    ta_m=0.5*(np.max(ta)+np.min(ta))
    hr=np.mean(hr)
    ea=es(ta)*(hr/100)
    ea_m=0.5*(np.max(ea)+np.min(ea))
    vpd_m=es(ta_m)-ea_m
    
    
    fhs=0.9*fhs+0.1
    fed=0.34-0.14*np.sqrt(ea_m)
    ft=4.9e-9*(ta_m+273.15)**4
    
    rb=fhs*fed*ft# %(MJ m-2 day-1)long wave radiation
    rs=(1-0.23)*Rs_day#(MJ m-2 day-1)short wave radiation
    rn=rs-rb
    
    delta=4098*es(ta_m)/(ta_m+237.3)**2#(kPa K-1)slope of the vapour pressur curve
    
    eto=((delta*rn+0.5*vpd_m*U2)/(2.45*(delta+0.067*(1+0.332*U2))))
    
    return eto
    
def es(Ta):
    """
    Compute saturated evaporation using Tetens formula (Campbell&Norman 1999, eq 3.8)
    
    Inputs
    -------
    Ta=air temeprature (C)
    
    Return
    -------
    es
    """
    a=0.611# constant (kPa)
    b=17.502#constant
    c=240.97#constant(C)
    
    es=a*np.exp(b*Ta/(Ta+c))
    
    return es


    