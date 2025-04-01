# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 06:07:23 2019

@author: Omar
"""

import numpy as np
import warnings

def tcanopy(rint,Ta,u_htop,Hr,Pa,h_top,Gc,G,psiH=0,psiM=0,
            epsilon_s=0.98):
    """
    Function to compute canopy temperature based on Campbell and Normam (1999)  
    Wind velocity is calculated at the top of the canopy.
    
    V1: Adiabatic profile is considered. Hence, no stability corrections are 
    considered. Tcanopy will be computed during the day. 
    During the night, air temperature will be used instead.
    
    Inputs:
    ---------
       rint= short wave net radiation (W/m2)
       Ta=air temperature (C)
       u_htop=wind velocity at the top of the canopy (m/s)
       Hr=relative humudity (%)
       Pa= atmospheric pressure (kPa)
       h_top=top of the canopy(m)
       Gc=Canopy conductance (mol m-2 s-1)
       G=soil sensible heat (w m-2)
       LAI=Leaf Area Index (m2 leaf m-2 soil)
       lw
       epsilon_s=Surface emissivity(dimensionless)
       boundary_layer= If false assumes a coupled canopy 
       psiM and psiH= stability correction if zero then adiabatic temperature 
                      profiles 
       
    Return:
    ---------
    Tcanopy= Canopy temperature (C)
       
    
    References:
    -------------
       Campbell G.S. and Norman J.M.1999. Plant and plant communities. Ed:
       Campbell G.S. and Norman J.M.1998.An introduction to environmental 
           Biophysics. pg:223-246
    
       Villalobos, F.J;Orgaz, F; Testi,L and Fereres,E. 2000. Measurement and
           modelling of evapotranspiration of olive orchards. EurJourAgr.pg155-163
           
    """
    #Constants
    cp=29.3#(J mol-1 K-1) specific heat of air
    cp_kJ=29.3e-3#(kJ mol-1 K-1) specific heat of air
    lambda_heat=44#(kJ mol-1)latent heat of vaporization
    Stef_Boltz=5.6697e-8#( W m-2 C-4) Stefan Boltzman constant
    gamma=cp_kJ/lambda_heat#(C-1) Thermodynamic psychrometer constant
    
    
    
    #Computing saturation slope
    b=17.502#Constant to compute saturation slope
    c=240.97#Constant to compute saturation slope
    es_t=es(Ta)
    s=saturation_slope(Ta,b,c,es_t,Pa)
    
    vpd=es(Ta)*(1-Hr/100)
    
    #Radiative conductance
    gr=(4*epsilon_s*Stef_Boltz*(Ta+273.15)**3)/cp
    
    #Heat and vapour boundary layer conductances. (The equation is the same for
    #both)
    
    rho_den=44.6*Pa/101.3*273.15/(Ta+273.15)#(mol m-3) gas molar density variation 
                                            #with temperature       
        
    d=0.65*h_top#zero plane displacement. Relation obtained for maize (Chapter 5, pg71)
    zm=0.15*h_top#momentum roughness parameter. The coefficient has been adjusted 
                 #based on estimations by Villalobos et al(2000)
    zh=0.2*zm

                    
   
    gHa=(0.4**2*rho_den*u_htop)/((np.log((h_top-d)/zm)+psiM)*(np.log((h_top-d)/zh)+psiH)) #Boundary layer conductance for heat(mol m-2 s-1) above field
    gva=gHa#Boundary layer conductance for vapour (mol m-2 s-1)
    gHr=gHa+gr#Sum of boundary layer and radiative conductance

    gv=1/(1/gva+1/Gc)#Conductance for vapour including boundary layer conductance (mol m-2 s-1)
 
    gamma_star=gamma*gHr/gv#Apparent psychrometer constant
    
            
    A=(rint-epsilon_s*Stef_Boltz*Ta**4-G)/(cp*gHr)
    
    B=vpd/(gamma_star*Pa)
    
    Tcanopy=Ta+gamma_star/(s+gamma_star)*(A-B)
        
    return Tcanopy
    
    
def saturation_slope(Ta,b,c,es,Pa):
    """
    Slope of the vapur saturation curve from Campbell and Norman (1999).
    Chapter 3.pg 41 
    
    Inputs:
    ---------
    Ta= air temperature (C)
    b=constant usually fixed at 17.502
    c= constant usually fixed at 240.97 (C)
    Pa= atmospheric pressure (kPa)
    
    Return
    --------
    s=saturation slope
    
    Reference:
    ----------
    Campbell, G. S. and J. M. Norman (1998). Introduction to environmental
        biophysics. New York, Springer.
    
    """

    delta=(b*c*es)/(c+Ta)**2

    s=delta/Pa
    
    return s

def es(Ta):
    """
    Compute saturated evaporation using Tetens formula (Campbell&Norman 1999, eq 3.8)
    
    Inputs
    -------
    Ta=air temeprature (C)
    
    Return
    -------
    es: Vapuro pressure at saturation (kPa)
    
    Reference:
    ----------
    Campbell, G. S. and J. M. Norman (1998). Introduction to environmental
        biophysics. New York, Springer.
    
    """
    a=0.611# constant (kPa)
    b=17.502#constant
    c=240.97#constant(C)
    
    es=a*np.exp(b*Ta/(Ta+c))
    
    return es

def G_surface(Ta,kc,kd,t,t0=8):
    """
    Function to calculate heat flux density at the soil surface.
    Calculations are based on Campbell and Norman (1999)
    
    Inputs:
    --------
       Ta=Air temperature (ºK)
       kc=Soil thermal conductivity(W m-1 K-1)
       kd=soil thermal diffusivity (m2 s-1)
       t=time of the day (s)
       t0=8;%Phase shift. It is set to 8 if local time is used.
       
    Return:
    -------
        G= Heat flux density at the soil surface (W m-2)
    
    References:
    -----------    
       Campbell G.S. and Norman J.M. Heat Flow in the soil.1999.Chap8 Ed:
       Campbell G.S. and Norman J.M  An introduction to environmental Biophysics. 
    """ 
    t0=t0*3600
    A_0=((max(Ta)-min(Ta))/2)#Amplitude of temperature fluctuations
    
    omega_d=2*np.pi/(24*3600)#Angular frequency daily basis
    
    D=np.sqrt(2*kd/omega_d)#Damping depth
    
    G=(np.sqrt(2)*kc*A_0*np.sin(omega_d*(t-t0)+np.pi/4))/D
    
    return G

def u_topcanopy(u,aneHeight,h_top):
    
    """
    Wind velocity at the top of the canopy
    
    Inputs:
    --------
       u: wind velocity at the anemometer (m s-1)
       aneHeight: height of the anemometer (m)
       h_top: canopy height (m)
       
    Return:
    -------
        u_htop= wind velocity at the top of the canopy (m s-1)
    
    References:
    -----------    
       Villalobos, F. J. and Elias F.C. (2017). 
       Fitotecnia: principios de agronomía para una agricultura sostenible. 
       Spain, Ed. Mundi-Prensa.(Eq 4.10)
    """ 
    

    z=h_top#height of the anemometer (m)
    if z!=2.0:
       #find velocity at 2m height using the reference u
       hmeadow=0.12#height of the reference medow
       dmedow=0.65*hmeadow#zero plane displacement. Relation obtained for maize (Chapter 5, pg71)
       zmedow=0.15*hmeadow#momentum roughness parameter. The coefficient has been adjusted based on estimations by Villalobos et al(2000)
       
       u2=u*(np.log(2-dmedow)-np.log(zmedow))/(np.log(z-dmedow)-np.log(zmedow))
        
    elif z==2.0:
        u2=u
        

    u_htop=2.6*u2/(6.6-np.log(h_top))
    
    return u_htop
    
    
def stab_corr(u,T,Tc,h_top,Pa,aneHeight):
    
    """
    Calculation of the stability correction factors for boundary layer 
    conductance determination
    
    Inputs
    ------
    u= wind velocity (m/s)
    T= air temperature (C)
    Tc= canopy temperature (C)
    h_top= canopy height (m)
    h_base= canopy base. From soil to the first green branch (m)
    lw= leaf width (m)
    LAI=leaf area index (m2leaf m-2soil)
    Pa= atmospheric pressure (kPa)
    aneHeight=anemometer height (m)
    
    Output
    ------
    psiM
    psiH
    
    Reference
    ---------
    Bitelli, M; Campbell, G.S.; Fausto. T.2015.Atmospheric Boundary conditions.
    Soil physiscs with basics (Ch 15).Bitelli, M; Campbell, G.S.; Fausto. T. 
    (Ed) Oxford University Press. pp:367-378
    """
    
    
    
    g=9.81 #gravity acceleration (m s-2)
    cp=29.3#(J mol-1 K-1) specific heat of air
    vk= 0.4 # von Kármán's constant 
    
    Tk=Tc+273.15# Surface temperature in Kelvins
    ro=44.6*(Pa/101.3)*(293.15/Tk)#molar density of the gas
    Ch=ro*cp #volumetric heat of air (=1200 J m-3 K-1 at 20 C and at sea level) 
    
    d=0.65*h_top#zero plane displacement. Relation obtained for maize (Chapter 5, pg71)
    zm=0.15*h_top#momentum roughness parameter. The coefficient has been adjusted based on estimations by Villalobos et al(2000)
    zh=0.2*zm
    
    z=h_top
    
    # if z<=h_top:
    #     #Find velocity at the top canopy. Upscaling is computed assumed that
    #     #anemometer is placed in a bare soil
    #     z_ane=h_top+2#Virtual anemometer height at 2m above the canopy (m)
    #     zm_s=0.002#Rougness at soil surface(m). Form Campbell&Norman(1999), table5.1
    #     u_star_ane=0.4*u/np.log(z/zm_s)#friction velocity at anemometer height (m s-1)
    #     u_ane=u_star_ane/0.4*np.log((z_ane)/zm_s)#Wind velocity corrected at a height equal to h_top
        
    #     u_star=0.4*u_ane/np.log((z_ane-d)/zm)#u_star with corrected velocity(m s-1)
    #     u_htop=u_star/0.4*np.log((h_top-d)/zm)#wind velocity at the top of the canopy (m s-1)
       
    # else:
    #     u_star=0.4*u/np.log((z-d)/zm)#friction velocity(m s-1)    
    #     u_htop=(u_star/0.4)*np.log((h_top-d)/zm)#Wind velocity at the top of the canopy(m/s)

    # h_canopy=(h_top+h_base)/2#Height at the canopy center (m)
    # h_canopy=zm+d#Height at the heat source/sink height
    # lm=np.sqrt(4*lw*h_top/(np.pi*LAI))#Mid distance between leaves
    # a=np.sqrt((0.2*LAI*h_top)/lm)#attenuation coefficient
    u_hmid=u#*np.exp(a*(h_canopy/h_top-1))#wind speed is exponentially attenuated from the top of the canopy to height of interest "h_canopy"

    
    psiM=0
    psiH=0
    psiH_positive=0
    psiH_negative=0
    errorLim=0.01
    error=50
    iteration=0
    while error >= errorLim:
                
        u_star=vk*u_hmid/(np.log((z-d+zm)/zm)+psiM)
        K=vk*Ch*u_star/(np.log((z-d+zh)/zm)+psiH)
        H=K*(T-Tc)#sensible heat
        zeta=-vk*z*g*H/(Ch*Tk*np.power(u_star,3))
        
        if zeta>0:
            psiHold_positive=psiH_positive
            psiH_positive=6*np.log(1+zeta)
            if psiH_positive!=0:
                error=np.abs((psiH_positive-psiHold_positive)/psiH_positive)*100
            psiH=psiH_positive
            psiM=psiH
        else:
            psiHold_negative=psiH_negative
            psiH_negative=-2*np.log((1+np.sqrt(1-16*zeta))/2)
            if psiH_negative!=0: 
                error=np.abs((psiH_negative-psiHold_negative)/psiH_negative)*100
            psiH=psiH_negative
            psiM=0.6*psiH
            
        if iteration>100:
            break
        else:
            iteration=iteration+1
        
    return psiM,psiH


def loop_stab(rint,Ta,u,Hr,aneHeight,Pa,h_top,gs,G,
                            epsilon_s=0.98):
    
        TcLower=Ta-15 #Lower limit for loop
        TcUpper=Ta+15 #Upper limit for loop
        Tcanopy=TcLower
        errorlim=0.01# error limit to stop iteration
        iteration=0
        error=50
        while error>=errorlim:   
            Tcanopy_old=Tcanopy              
            Tcanopy=(TcLower+TcUpper)/2
            
            
            
            
            if Tcanopy != 0:
                error=np.abs((Tcanopy-Tcanopy_old)/Tcanopy)*100
                
            
            psiMlower,psiHlower=stab_corr(u,Ta,TcLower,
                                          h_top,Pa,aneHeight)
            
            TcLower=tcanopy(rint,Ta,u,Hr,Pa,h_top,
                            gs,G,psiH=psiHlower,psiM=psiMlower,
                            epsilon_s=epsilon_s)
            
            
            psiMupper,psiHupper=stab_corr(u,Ta,TcLower,
                                          h_top,Pa,aneHeight)
            
            TcUpper=tcanopy(rint,Ta,u,Hr,Pa,h_top,
                            gs,G,psiH=psiHupper,psiM=psiMupper,
                            epsilon_s=epsilon_s)
            
            test=TcUpper*TcLower
            
            if test<0:
                TcUpper=Tcanopy
            elif test>0:
                TcLower=Tcanopy
            else:
                error=0                    
                
                
            if iteration>100:
                Tcanopy=tcanopy(rint,Ta,u,Hr,aneHeight,Pa,h_top,
                            gs,G,psiH=0,psiM=0,
                            epsilon_s=epsilon_s)
                break
            
            iteration=iteration+1
            
        return Tcanopy   
            
def tc_stab (rint,Ta,u,Hr,Pa,h_top,Gc,G,epsilon_s=0.98): 


      g=9.81 #gravity acceleration (m s-2)
      cp=29.3#(J mol-1 K-1) specific heat of air
      vk= 0.4 # von Kármán's constant 
    
      
      d=0.65*h_top#zero plane displacement. Relation obtained for maize (Chapter 5, pg71)
      zm=0.15*h_top#momentum roughness parameter. The coefficient has been adjusted based on estimations by Villalobos et al(2000)
      zh=0.2*zm
    
      z=h_top
      
      #Initial conditions
      psiM=0
      psiH=0
      Tc=Ta
      error=100
      errorlim=0.01
      niteration=0
      while error>=errorlim:
          Tc_old=Tc
          psiM_old=psiM
          psiH_old=psiH
          Tk=Tc_old+273.15# Surface temperature in Kelvins
          ro=44.6*(Pa/101.3)*(293.15/Tk)#molar density of the gas
          Ch=ro*cp #volumetric heat of air (=1200 J m-3 K-1 at 20 C and at sea level) 
        
          u_star=vk*u/(np.log((z-d+zm)/zm)+psiM_old)
          K=vk*Ch*u_star/(np.log((z-d+zh)/zm)+psiH_old)
          H=K*(Ta-Tc_old)#sensible heat
          zeta=-vk*z*g*H/(Ch*Tk*np.power(u_star,3))
          if zeta>0.0:
            psiH_positive=6*np.log(1+zeta)  
            psiH=psiH_positive
            psiM=psiH
          elif zeta<0.0:
            psiH_negative=-2*np.log((1+np.sqrt(1-16*zeta))/2)
               
            psiH=psiH_negative
            psiM=0.6*psiH
            
          elif zeta==0.0:
              psiH=0
              psiM=0
           
          Tc=tcanopy(rint,Ta,u,Hr,Pa,h_top,Gc,G,psiH=psiH,psiM=psiM,
                epsilon_s=0.98) 
          
          error=np.abs((Tc-Tc_old)/Tc)*100
          
          if niteration>100:
              break
          niteration=niteration+1
          
      return Tc
         
            
            
            
    
# if __name__=='__main__':
    
#         u_htop=np.linspace(0.8,10)
        
# #       #Tc2=tcanopy(400,40,u_htop,20,101,2,0.4,100,psiH=0,psiM=0,
# #       #              epsilon_s=0.98)
      
#         Tc=np.zeros(u_htop.shape[0])
#         i=0
#         for u in u_htop:
#             Tc[i]=tc_stab(400,40,u,20,101,2,0.12,100,
#                               epsilon_s=0.98)
          
#             i=i+1
    
#     gs=np.linspace(0.001,0.4)
#     Tc=np.zeros(gs.shape[0])
#     i=0
#     for g in gs:
        
#         Tc[i]=tcanopy(300,31.8,4.1,28,2,98.85,2.5,g,67.81,psiH=0,psiM=0,
#                 epsilon_s=0.98, boundary_layer=True)
#         i=i+1
#      stab_corr(2,35,40,3,1,0.08,6,101,2)