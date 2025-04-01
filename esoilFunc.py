# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 23:36:10 2020

@author: Omar
"""

import numpy as np

def delta(ta):
    
    """
    Slope of the saturation curve
    
    Input:
    ------    
        ta: air temperature (C)
        
    Output:
    -------    
        delta: slope of the vapor pressure curve(kPa C-1)
        
    Reference:
    ----------    
        Allen,R.G.,Pereira,L.S.,Raes,D.and Smith, M.1998.
        Crop evapotranspiration- Guidelines for computing crop water 
        requirements. In "FAO irrigation and Drainage paper", 
        Vol. 56, pp. 322. 
        Food and Agriculture Organization of the United Nations, Rome.
    """
    
    delta= 4098*es(ta)/(ta+237.3)**2
    
    return delta

def gamma(patm):
    
    """
    Psicrometric constant
    
    Input:
    ------    
        patm: atmospheric pressure
        
    Output:
    -------    
        gamma: psicrometric constant (kPa C-1)
        
    Reference:
    ----------    
         Allen,R.G.,Pereira,L.S.,Raes,D.and Smith, M.1998.
        Crop evapotranspiration- Guidelines for computing crop water 
        requirements. In "FAO irrigation and Drainage paper", 
        Vol. 56, pp. 322. 
        Food and Agriculture Organization of the United Nations, Rome.
        
    """
    
    gamma=0.000665*patm
    
    return gamma

def patm(alt):
    
    """
    Atmospheric pressure variation with altitude
    
    Input:
    ------    
        alt: altitude (m.a.s.l)
        
    Output:
    -------    
        patm: atmospheric pressure at the current altitude (kPa)
        
    Reference:
    ----------    
        Villalobos, F. J., and Elias, F. C. (2017). 
        "Fitotecnia: principios de agronomía para una agricultura sostenible," 
        Ed. Mundi-Prensa, Spain.
    """
    
    patm=101.3*(1-alt/44308)**5.2568
    
    return patm

def eso(delta,gamma,rnsoil,vpd,u,fwet):
    
    """
    Potential soil evaporation at the energy limiting stage in the area
    wetted by the emitter
    
    Input:
    ------    
        delta: solpe of the saturation vapor pressure (kPa C-1)
        gamma: psicrometric constant (kPa C-1)
        rnsoil: net radiation reaching the soil (MJ m-2 day-1)
        vpd: average daily vapor pressure deficit (kPa)
        u: average daily wind velocity (m s-1)
        fwet: fraction of wetted area by the emitter (dimensionless)
        
    Output:
    -------    
        eso: soil evaporation during the limiting stage (mm day-1)
        
    Reference:
    ---------    
        Bonachela, S., Orgaz, F., Villalobos, F., and Fereres, E. 
        (1999). Measurement and simulation of evaporation from soil in olive 
        orchards. Irrigation Science 18.
        
        Allen,R.G.,Pereira,L.S.,Raes,D.and Smith, M.1998.
        Crop evapotranspiration- Guidelines for computing crop water 
        requirements. In "FAO irrigation and Drainage paper", 
        Vol. 56, pp. 322. 
        Food and Agriculture Organization of the United Nations, Rome.
    """
    
    rnsoil=0.408*rnsoil#Net radiation in mm day-1 (FAO 56)
    #microadvective coefficient (dimensionless)
    if fwet>0:
        rnsoil=rnsoil*fwet
        kw=1.1+0.14*np.log(1/fwet)
    else:
        kw=1
    eso= (delta/(delta+gamma)*rnsoil+gamma*(delta+gamma)*vpd*2.7*(1+u/100))*kw
    
    return eso

def Esoil(ETo,ta,patm,Qd,rnsoil,vpd,u,fwet,fdry,CEwini,CEdini,
          twet,tdry,Ue=8,ce=3.5):
    
       """
        Soil evaporation from dry and wetted areas based on Ritchie. Potential 
        evaporation from the wetted area is based on Bonachela (2001). 
        
        Inputs:
            
            ETo: reference evapotranspiration
            ta: air temperature (C)
            patm: atmospheric pressure (kPa)
            Qd: Fraction of intercepted radiation
            rnsoil: net radiation reaching the soil (W m-2)
            vpd: vapor pressure deficit (kPa)
            u: wind velocity (m s-1)
            fwet: fraction of wetted area by the emitter (dimensionless)
            fdry: fraction of dry area (dimensionless)
            CEwini:cummulative evaporation from the soil wetted area (mm)
            CEdini:cummulative evaporation from the soil dry area (mm)
            twet:time since last wetted event (irrigation or precipitation) in the
                 wetted area (day)
            tdry:time since last wetted event (precipitation) in the
                 dry area (day) 
            Ue: coefficient for transition from phaseI to phaseII (mm)
            ce:shape parameter for phaseII evaopration (mm day-0.5)
            
        Output:
            
            Eswet: evaporation from the soil wetted area (mm day-1)
            Esdry: evaporation from the soil dry area (mm day-1)
            CEwet: updated cummulative evaporation from the wetted area
            CEdry: updated cummulative evaporation from the dry area
                 
        
        
        
        
       """
       
       
       DELTA=delta(ta)
       GAMMA=gamma(patm)
       if CEwini <= Ue and CEdini>Ue:           
           Eswet=eso(DELTA,GAMMA,rnsoil,vpd,u,fwet)
           
       elif CEwini<=Ue and CEdini<=Ue:
           Eswet=ETo*(1-Qd)*fwet
           
       else:
           Eswet=ce*(twet**0.5-(twet-1)**0.5)
           
       if CEdini<=Ue:
           Esdry=ETo*(1-Qd)*fdry
       else:
           Esdry=ce*(tdry**0.5-(tdry-1)**0.5)
        
       CEwet=CEwini+Eswet
       CEdry=CEdini+Esdry
       
       return Eswet,Esdry,CEwet,CEdry
   
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