# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 06:59:41 2019

@author: Omar
"""

"""
SPAC model functions

This class contains the external functions using for SPAC model calculations 
except the ones related to the radiation calculations
 
"""


import numpy as np
from scipy.optimize import newton
import warnings
##############################PHOTOSYNTHESIS FUNCTIONS########################

def tfunc(ta,cpar,delh):
    """
     Calculation of the prarameters of Farquar equation under temperature
     limitation (Bernachi,2001)
     
     Input:
     ------    
         ta: air temperature (C)
         cpar: scaling constant
         delh:activation energy
         
     Output:
     -------   
         y= adjusted farquhar parameter to temperature
    
    Reference:
        
        Bernacchi, C.J.; Singsaas, E.L.; Pimentel, C.; Portis Jr, A.R.; 
        Long, S.P. 2001. Improved temperature response functions for models of
        Rubisco-limited photosynthesis. Plant, Cell & Environment. pp 253-259
     
     
    """
    tak=ta + 273 #Ambient temperature in Kelvins
    R=8.31 #gas constant (J/K mol)
    y = np.exp(cpar-delh/(R*tak))
    return y

def solve_tuzet_joint_energy(ca,zalfa1,zalfa2,LAI,tleaf,ta,Hr,gva,g0p,mp,sf,potr,
                      ab_t, bb_t, eb_t, db_t,ab_r, bb_r, eb_r, db_r, tes,rd,Pa):
     """
    Iterative approach to find equilibrium solution for stomatal conductance,
    net assimilation, internal CO2 concentrarion and leaf water potential.This 
    apprach was developed by Prof Francisco J. Villalbos during the SPAC model 
    development.
    
    In this "energy version" the boundary layer of the leaf is considered. 
    
    Inputs:
    -------   
        ca: ambient CO2 concentration (micromol CO2 mol air-1)
        zalfa1: resistance of the soil
        zalfa2 resistance of teh canopy
        LAI:leaf area index (m2 leaf m2 soil-1)
        tleaf: leaf temeprature (C)
        ta: Air temperature (C)
        Hr:relative Humidity
        gva: Leaf boundary layer conductance for water vapour (mol m-2 s-1)
        g0p: Residual stomatal conducntance at compensation point (mol m-2 s-1)
        mp: Tuzet's adjusting parameter (dimensionless)
        potr: Leaf water potential at which stomatal conductance start to be 
              limited (kPa)
        sf: Shape parameter (-kPa)
        Pa:atmospheric pressure (kPa)
        tes: CO2 compensation point (micromol CO2 mol-1) 
        rd: Rate of CO2 evolution in the light resulting from process other
            than photorespiration (micromol CO2 m-2 s-1) 
        
        * The rest of the inputs are algebraic representations of the equations
        to compute photosynthesis when Rubisco or RuBp is limited
        
    
     Outputs:
    -------    
        gs: stomatal conductance for CO2 (mol m-2 s-2)
        evap: leaf transpiration assuming no boundary layer conductance(gH2O m-2 s-1)
        asim: net assimilation (micromol m-2 s-1)
        psleaf: leaf water potential. Referred as the water potential at the 
                petiole (kPa)
        ci: internal CO2 concentration (micromol CO2 mol-1 air)
        
        * All the m2 are referred to leaf area
    
    References:
    -----------
        Farquhar, G.D.; von Caemmerer, S.; Berry, J.A. 1980. A biochemical model
       of photosynthetic CO2 assimilation in leaves of C3 species. Planta pp 78-90
       
       Tuzet, A.; Perrier, A.; Leuning, R. 2003. A coupled model of stomatal 
       conductance, photosynthesis and transpiration. Plant, Cell & Environmnent.
       pp 1097-1116
       
       Campbell, G.S. and Norman J.M. 1998. An introduction to environmental
            biophysics. Springer. New York. pp 286.Chapter 14.    
        
     """
     #Equilibrium ci  
     ci_eq=(ab_r* bb_t - ab_t * bb_r) / (ab_t * eb_r - ab_r * eb_t)
     
     #Atmosphere vapour pressure 
     ea=es(ta)*Hr/100
     
     ci_sol=999
     
     if np.isnan(ca):
         raise Exception('Not CO2 value. Please check weather file')
         
     for ci in np.arange(0.98*ca,0.1*ca,-1.0):   
        if ci<ci_eq:
            gs = (ab_t*(ci-tes)-rd*(eb_t*ci+bb_t))/((eb_t*ci+bb_t)*(ca-ci))#conductance temp limited(mol/ m2/s)
        
        else:
            gs = (ab_r*(ci-tes)-rd*(eb_r*ci+bb_r))/((eb_r*ci+bb_r)*(ca-ci))#conductance rad limited(mol/ m2/s)
            
           
        
        if gs>0.8:
            gs=0.8
        
        
        asim = gs*(ca-ci) #net asimilation (micromol/m2/s)
        
        gsw=1.6*gs
        gv=gvap(gva,gsw)
        psileaf = zalfa1-(zalfa2)*gv*(es(tleaf)-ea)/Pa*0.018*LAI#the factor 18e-3 converts
                                                      #evaporation in mol of H2O to kg. 
                                                      #Units of psileaf (kPa)
    
        
    
        fwat = (1 + np.exp(sf*potr))/(1+np.exp(sf*(potr - psileaf)))
        gs2 = g0p+(mp*(asim+rd)/(ci-tes))*fwat                    
        
        dife=gs-gs2
        if ci<0.98*ca:    
                if np.sign(dife)!=np.sign(dife_prev):
                   ci_sol=ci-dife*(ci_prev-ci)/(dife_prev-dife)
                   break
                
        
        ci_prev=ci
        dife_prev=dife
       
     
       
     ci=ci_sol
            
     if ci==999:       
         
          for ci in np.arange(1.01*ca,5*ca,1.0):   
            if ci<ci_eq:
                gs = (ab_t*(ci-tes)-rd*(eb_t*ci+bb_t))/((eb_t*ci+bb_t)*(ca-ci))# conductance temp limited(mol/ m2/s)
    
            else:
                gs = (ab_r*(ci-tes)-rd*(eb_r*ci+bb_r))/((eb_r*ci+bb_r)*(ca-ci))#conductance rad limited(mol/ m2/s)
    
            
            asim = gs*(ca-ci)#net asimilation (micromol/m2/s)
            
            gsw=1.6*gs
            gv=gvap(gva,gsw)
            psileaf = zalfa1-zalfa2*gv*(es(tleaf)-ea)/Pa*0.018*LAI#the factor 18e-3 converts evaporation in mol of H2O to kg. Units of psileaf (kPa)
            
            fwat = (1 + np.exp(sf*potr))/(1+np.exp(sf*(potr - psileaf)))          
            gs2 = g0p+(mp*(asim+rd)/(ci-tes))*fwat                    
    
            dife=gs-gs2
            if ci>1.02*ca:    
                if np.sign(dife)!=np.sign(dife_prev):
                   ci_sol=ci-dife*(ci_prev-ci)/(dife_prev-dife)
                   break
                
            
            ci_prev=ci
            dife_prev=dife
          
           
          
         
      
     ci=ci_sol       
     if ci<ci_eq:
        gs = (ab_t*(ci-tes)-rd*(eb_t*ci+bb_t))/((eb_t*ci+bb_t)*(ca-ci))#conductance temp limited(mol/ m2/s)
    
     else:
         gs = (ab_r*(ci-tes)-rd*(eb_r*ci+bb_r))/((eb_r*ci+bb_r)*(ca-ci))# conductance rad limited(mol/ m2/s)
    
       
        
     if gs<g0p:
        gs=g0p
        asim=(ca-ci)*gs
        
        gsw=1.6*gs
        gv=gvap(gva,gsw)
        psileaf = zalfa1-(zalfa2)*gv*(es(tleaf)-ea)/Pa*0.018*LAI#the factor 18e-3 converts evaporation in mol of H2O to kg. Units of psileaf (kPa)
     
     return ci,asim,gs,psileaf
 
def tuzet_energy(ta,Rint,ca,u,Hr,alt,leaf_w,
                 zalfa1,zalfa2,LAI,g0p,mp,sf,potr,
                 THETA,ALFA,CTES,DELHTES,ckc,DELHKC,cKO,
                 DELHKO,CRD,DELHRD,CVCX,DELHVCX,CJX,
                 DELHJX,CTPU,DELHTPU):
    
    """
    Solve transpiration, condunctance, assimilation and Tleaf applying the 
    energy balance clousure.
    The calculations of Tleaf are based on the linearization proposed by 
    Campbell and Norman (1998) Eq.(14.6)
    
    
    
    """
    #Constant
    cp=29.3#(J mol-1 K-1) specific heat of air
    Mwg=18#molecular mass of water (g/mol)
    Pa = 101.3 * (1 - alt / 44308) ** 5.2568
    lambda_heat=44e3#(J mol-1)latent heat of vaporization
    ea=es(ta)*Hr/100
    
    
    #Convert PAR Radiation to microeinstein/m2/s assuming an averge
    #wave length fro PAR radiation of 550nm (Campbell&Norman, 1998)
    
    rint_par=Rint/0.217
    
    
    OI=210#intercellular O2 fraction
    
    
    
    gva=gvabl (u,leaf_w)
    gHa=gHabl(u,leaf_w)
    tl=ta
    error_tole=0.001
    error=100
    niter=0
    while error>error_tole:
        tl_old=tl
        
    
        #Farquhar parameters corrected for temperature
       
        tes = tfunc(tl_old,CTES, DELHTES)#CO2 compensation point in the absence of Rd (no units)
        KC = tfunc(tl_old,ckc, DELHKC)#Michaelis constant for CO2 (no units)
        KO = tfunc(tl_old,cKO, DELHKO)#Michaelis constant for O2 (no units)
        rd = tfunc(tl_old,CRD, DELHRD)#Respiration (no units)
        VCMAX = tfunc(tl_old,CVCX, DELHVCX)#maximum catalytic activity of Rubisco in the presence of saturating amounts of RuBP and CO2 (no units)
        JMAX = tfunc(tl_old,CJX, DELHJX)#maximum rate of electron transport at saturating irradiance (no units)
    #    TPU = TFUNC(CTPU, DELHTPU, tak)#rate of Pi release associated with triose phosphate utilization ()
        
        #parameters for temperature limitation
        ab_t = VCMAX
        db_t = tes
        eb_t = 1
        bb_t = KC*(1+OI/KO)
        
        #parameters for rad limitation
        a = THETA
        b = -(ALFA * rint_par + JMAX)
        c = ALFA * rint_par * JMAX
        j = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
        
        ab_r = j
        db_r = tes
        eb_r = 4
        bb_r = 8 * tes
        
        ci,asim,gs,psileaf=solve_tuzet_joint_energy(ca,zalfa1,zalfa2,LAI,tl_old,
                                                    ta,Hr,gva,g0p,mp,sf,potr,
                                                    ab_t, bb_t, eb_t, db_t,ab_r, 
                                                    bb_r, eb_r, db_r, tes,rd,Pa)
    
        
        gsw=1.6*gs
        gv=gvap(gva,gsw)
       
        tl=Tleaf(Rint,u,ta,Hr,alt,leaf_w,gsw,gs_ad=0,epsilon_s=0.98)

        
        closure=lambda tleaf: Rint-cp*gHa*(tleaf-ta)-lambda_heat*gv*(es(tleaf)-ea)/Pa
        tl=newton(closure,tl)
        
        error=np.abs((tl-tl_old)/tl)*100
        
        
        niter=niter+1
        
        if niter>100:
            warnings.warn('Number of iterations exceeded')
            tes = tfunc(ta,CTES, DELHTES)#CO2 compensation point in the absence of Rd (no units)
            KC = tfunc(ta,ckc, DELHKC)#Michaelis constant for CO2 (no units)
            KO = tfunc(ta,cKO, DELHKO)#Michaelis constant for O2 (no units)
            rd = tfunc(ta,CRD, DELHRD)#Respiration (no units)
            VCMAX = tfunc(ta,CVCX, DELHVCX)#maximum catalytic activity of Rubisco in the presence of saturating amounts of RuBP and CO2 (no units)
            JMAX = tfunc(ta,CJX, DELHJX)#maximum rate of electron transport at saturating irradiance (no units)
        #    TPU = TFUNC(CTPU, DELHTPU, tak)#rate of Pi release associated with triose phosphate utilization ()
            
            #parameters for temperature limitation
            ab_t = VCMAX
            db_t = tes
            eb_t = 1
            bb_t = KC*(1+OI/KO)
            
            #parameters for rad limitation
            a = THETA
            b = -(ALFA * rint_par + JMAX)
            c = ALFA * rint_par * JMAX
            j = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
            
            ab_r = j
            db_r = tes
            eb_r = 4
            bb_r = 8 * tes
            
            ci,asim,gs,psileaf=solve_tuzet_joint_energy(ca,zalfa1,zalfa2,LAI,tl_old,
                                                        ta,Hr,gva,g0p,mp,sf,potr,
                                                        ab_t, bb_t, eb_t, db_t,ab_r, 
                                                        bb_r, eb_r, db_r, tes,rd,Pa)
                
            
            
            evap=1.6*gs*(es(ta)-ea)/Pa*Mwg
            break
        else:
    
            evap=gv*(es(tl)-ea)/Pa*Mwg
            
            
    
    return gs,asim,psileaf,evap,ci,tl
    
    
    
def leuning (g0p,mp,D0,ci,asim,tau,vpd):
    """
    

    Parameters
    ----------
    g0p : float
         Residual stomatal conducntance at compensation point (mol m-2 s-1)
        
    mp : float
        Convert assimilation to gs units (dimensionless).
    D0 : float
        Stomata sensitivity to VPD.
    ci : float
        Sub-stomata CO2 concentration.
    asim : float
        Net assimilation.
    tau : float
        CO2 compensation point.
    vpd: float
        vapor pressure deficit (kPa)

    Returns
    -------
    gs: float
        Stomata conductance for CO2 (mol m-2 s-1).
        
    Reference:
    ---------    
    Leuning, R. (1995). "A critical appraisal of a combined 
    stomatal-photosynthesis model for C3 plants." 
    Plant, Cell & Environment 18(4): 339-355.

    """    
    gs = g0p+mp*asim/((ci-tau)*(1-vpd/D0))
    
    return gs


    
    
def gHabl (u,leaf_w):
    
    """
    Boundary layer conductance for heat
    
    Input:
    ------    
        u: wind velocity (m s-1)
        leaf_w: leaf witdth (m)
               
    Output:
    -------    
        gHa: conductance for heat (mol m-2 s-1) 
                
    Reference:
    ----------    
        Campbell, G.S. and Norman J.M. 1998. An introduction to environmental
            biophysics. Springer. New York. pp 286.Chapter 14
    
    
    """
    d=0.7*leaf_w
    gHa=1.4*0.135*np.sqrt(u/d)
    
    return gHa

def gvabl (u,leaf_w):
    
    """
    Boundary layer conductance for water vapour
    
    Input:
    ------    
        u: wind velocity (m s-1)
        leaf_w: leaf witdth (m)
               
    Output:
    -------    
        gva: boundary layer conductance for water vapour (mol m-2 s-1) 
                
    Reference:
    ----------    
        Campbell, G.S. and Norman J.M. 1998. An introduction to environmental
            biophysics. Springer. New York. pp 286.Chapter 14
    
    
    """
    d=0.7*leaf_w
    gva=1.4*0.147*np.sqrt(u/d)
    
    return gva

def gvap (gva,gs_ab,gs_ad=0):
    
    """
    Boundary layer plus stomata conductance for water vapour
    
    Input:
    ------    
        gva: Boundary layer conductance for water vapuor (mol m-2 s-1)
        gs_ab: Stomata conductance abaxial (mol m-2 s-1)
        gs_ad: Stomata conductance adaxial (mol m-2 s-1)
               
    Output:
    -------    
        gv: Leaf conductance for water vapour (mol m-2 s-1) 
                
    Reference:
    ----------    
        Campbell, G.S. and Norman J.M. 1998. An introduction to environmental
            biophysics. Springer. New York. pp 286.Chapter 14
    
    
    """
    
    gv=0.5*gs_ab*gva/(gs_ab+gva)+0.5*gs_ad*gva/(gs_ad+gva)
    #gv=1/(1/gs_ab+1/gva)
    return gv

def grad (Ta,epsilon_s=0.98):
    
    """
    Radiative conductance
    
    Input:
    ------    
        Ta: air temperature (C)
        epsilon_s=Surface emissivity(dimensionless)
               
    Output:
    -------    
        gr: Radiative conductance (mol m-2 s-1) 
                
    Reference:
    ----------    
        Campbell, G.S. and Norman J.M. 1998. An introduction to environmental
            biophysics. Springer. New York. pp 286.Chapter 14
    
    
    """
    Stef_Boltz=5.6697e-8#( W m-2 C-4) Stefan Boltzman constant
    cp=29.3#(J mol-1 K-1) specific heat of air
    
    
    gr=(4*epsilon_s*Stef_Boltz*(Ta+273.15)**3)/cp
    
    return gr
    
    

def Tleaf(rint,u,Ta,Hr,alt,leaf_w,gs_ab,gs_ad=0,epsilon_s=0.98):
    """
    Leaf temperature calculation based on the equation proposed in 
    Campbell&Norman (1998). Chapter 14 Eq 14.6. 

    Parameters
    ----------
    rint : float
        Radiation intercepted by the leaf (W m-2).
    u : float
        wind velocity (m s-1).
    Ta : float
        Air temperature (ºC).
    Hr : float
        Relative humidity (%).
    alt : float
        Altitude (m.a.s.l).
    leaf_w : float
        Leaf width (m).
    gs_ab : float
        Stomata conductance for water vapour abaxial (mol m-2 s-1).
    gs_ad : float, optional
        Stomata conductance adaxial (mol m-2 s-1). The default is 0.
    epsilon_s : optional, optional
        Surface emissivity(dimensionless). The default is 0.98.

    Returns
    -------
    Tl : float
        Leaf temperature ºC.
        
    Reference
    --------
    Campbell, G.S. and Norman J.M. 1998. An introduction to environmental
            biophysics. Springer. New York. pp 286.Chapter 14
    """
    
    #Constants
    cp=29.3#(J mol-1 C-1) specific heat of air
    cp_kJ=29.3e-3#(kJ mol-1 C-1) specific heat of air
    lambda_heat=44#(kJ mol-1)latent heat of vaporization
    
    gamma=cp_kJ/lambda_heat#(C-1) Thermodynamic psychrometer constant
    
    Pa = 101.3 * (1 - alt / 44308) ** 5.2568
    #Computing saturation slope
    b=17.502#Constant to compute saturation slope
    c=240.97#Constant to compute saturation slope
    es_t=es(Ta)
    s=saturation_slope(Ta,b,c,es_t,Pa)
    
    vpd=es(Ta)*(1-Hr/100)
    
    #Radiative conductance
    gr=grad(Ta)
    
    #Conductance for Heat
    gHa=gHabl (u,leaf_w)
    
    #Boundary layer conductance for water vapour
    gva=gvabl(u,leaf_w)
    
    #Conductance for water vapour
    gv=gvap(gva,gs_ab)
    
    gHr=gHa+gr
    
    gamma_star=gamma*gHr/gv
    
    Tl=Ta+gamma_star/(s+gamma_star)*(rint/(gHr*cp)-vpd/(Pa*gamma_star))
    
    return Tl
    
    
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
    """

    delta=(b*c*es)/(c+Ta)**2

    s=delta/Pa
    
    return s
 
 
##########################PLANT RESISTANCE FUNCTIONS###########################
        
def k_xylem(pstem,Kmax=0.004,B=3000,C=3):
    """
     Function to compute xylem conductivity Sperry et al (2016) Equation.
    The model account for cavitation effects.
    
    Inputs
    ------
      Kmax=maximum hydraulic contance of the xylem (kg kPa-1 m-1 s-1)
      pstem= Pressure at the xylem (kPa)Absoulte values
      B= presure at wich k_xylem/Kmax=0.37
      C= Shape parameter 
    
    Default Parameter for Kmax and b are for olives have been obtained from Hacke et al (2014)
    C parameter was obtained from Sperry et al (2016)
    
    Output
    ------
     k_xyl: xylem conductivity (kg kPa-1 m-1 s-1)
    
    References
    ----------
    
    """

    k_xyl=Kmax*np.exp(-(np.abs(pstem)/B)**C)
    return k_xyl





def Psoil(sw,swsat,pe,b):
    
        """
        Soil water potential
        
        Inputs:
        -------    
            sw: soil water content (m3 m-3)
            swsat:saturated soil water content (m3 m-3)
            pe: air entry water potetial (kPa)
            b: shape parameter
            
        Output:
        -------    
            psoil: soil water potential (kPa)
            
        References:
        -----------    
            Campbell, G.S. and Norman J.M. 1998. An introduction to environmental
            biophysics. Springer. New York. pp 286
            
        """
        
        psoil= pe*(sw/swsat)**-b
        
        
        return psoil

def Ksoil(pe,b,psoil,ks):
    """
        Soil water potential
        
        Inputs:
        -------    
            psoil: soil water potential (kPa)
            ks:saturated soil water conductivity (kg s m-3)
            pe: air entry water potetial (kPa)
            b: shape parameter
            
        Output:
        -------    
            ksoil: soil water conductivity (kg s m-3)
            
        References:
        -----------    
            Campbell, G.S. and Norman J.M. 1998. An introduction to environmental
            biophysics. Springer. New York. pp 286
            
    """
    
    kSoil=ks*(pe/psoil)**(2+3/b)
    
    return kSoil





def Rsoil(lv,aroot,kSoil,dz):
    
    """
    Resistance from the soil to the root surface
    
    Inputs:
    ------    
        lv: root length density (m root m-3 soil)
        aroot: root radious (m)
        k:  soil water conductivity (kg s m-3)
        dz: layer depth
        
    Outputs:
    --------
        rsoil: soil resistance (s m-2 kPa kg-1)
        
    Reference:
    ----------
     Gardner, W.R. 1960. Dynamic aspects of water availability to plants. Soil
     science. pp 63-73
     
    """
        
    rsoil=np.log(1/(np.pi*lv*aroot**2))/(4*np.pi*kSoil*lv*dz)
    
    
    return rsoil

def RrootT (ta,aroot, y0=64934,a=1.09e9,b=-3.64):
    
    """"
    Root radial resistance as a function of soil temperature.
    Default Parameters correspond to olive tree roots 
    
    Input:
    ------    
        ta: air temperature (C)
        aroot: root radious (m)
        y0: adjust parameter
        a: adjust parameter
        b: adjust parameter
        
    Output:
    -------    
        rrT: root radial resistance as a function of temperature ( kPa m s kg-1)
        
    Reference:
    ----------    
        Garcia-Tejera, O.;Lopez-Bernal, A.; Villalobos, Francisco J.;Orgaz, F.;
        Testi, L. 2016. Effect of soil temperature on root resistance: 
            implications for different trees under Mediterranean conditions.
            Tree physiology. pp 469-478
        
   
    
    """
    
    rrT=(y0+a*ta**b)/(2*np.pi*aroot)*1000
    
    return rrT

def RrootTheta(rrT,sw,swsat,alfa=2,beta=30,sigma=0.25):
    
    """
    Root radial resistance as a function of soil water content.
    Default Parameters are those for a 
    
    Input:
    ------    
        rrT:root radial resistance as a function of temperature ( kPa m s kg-1)
        sw: soil water content (m3 m-3)
        swsat:saturated soil water content (m3 m-3)
        alfa:root resistance value at which sw/ swsat = sigma (Dimensionless)
        beta:Rate at which root resistance approaches to infinity (Dimensionless)
        sigma:Critical value of sw/swsat at which root resistance becomes 
               limiting	(Dimensionless)
               
    Output:
    -------    
        rrTheta: Radial root resistance as a function of soil water content 
                (kPa m s kg-1)
                
    Reference:
    ----------    
        BristowBristow, K. L., Campbell, G. S., and Calissendorff, C.1984.
        The effects of texture on the resistance to water-movement within 
        the rhizosphere. Soil Science Society of America Journal 48, 266-270.
        
    """

    
    rrTheta=rrT*(1+alfa*np.exp(-beta*(sw/swsat-sigma)))
    
    return rrTheta

#######################WATER BALANCE FUNCTIONS#################################

def rainfall_part(LAI,gc,rain):
    """
    Computes precipitation intercepted by the canopy
     Function to compute the precipitation reaching the soil
     through the calculation of the intercepted by the canopy.
     
    Inputs:
    -------
        gc: Ground cover
        rain: Precipitation (mm)
        LAI: leaf area index (m2 leaf m-2 soil)
        
    Output:
    ------
        pre_soil: amount of rain reaching the soil (mm)
    """

    XLLAI = LAI / gc
        
    STOR = 0.49 * XLLAI + 1.2
    CFAC = 2.41 + 0.55 * XLLAI
    PREI = STOR * (1 - np.exp(-rain / CFAC)) * gc
         
    pre_soil=rain-PREI
    return pre_soil

def theta_residual(pe,b,thetaSat,psoil_residual=-1e6):
    """
    Residual soil water content
    
     Inputs:
        -------    
            
            ks:saturated soil water conductivity (kg s m-3)
            pe: air entry water potetial (kPa)
            b: shape parameter
            thetaSat: soil water content at saturation (m3 m-3)
            psoil: oven-dry soil matric potential (kPa)
            
        Output:
        -------    
            theta_res: residual soil water content (m3 m-3)
            
        References:
        -----------    
            Campbell, G.S. and Norman J.M. 1998. An introduction to environmental
            biophysics. Springer. New York. pp 286
            
            Bittelli, M., Campbell, G., and Tomei, F. (2015). 
            "Soil Physics with Python. 
            Transport in the Soil–Plant–Atmosphere System," 
            Oxford University Press, New York, United States of America. pp 486
            
    """
    
    theta_res=thetaSat*(pe/psoil_residual)**(1/b)
    
    return theta_res

def varLp(psoil,m,Lpmin,psoil50=-733.6,Lpmax=False):
    """
    Function to vary Lp according to soil water potential. The variation
    follows a logistic function. The function type has been selected according
    to the results reported in Caldeira et al (2014).

    Parameters
    ----------
    psoil : float
        Soil water potential (kPa).
    m : float
        Rate of change of Lp.
    Lpmin : float
        Minimum root radial conductance (Kg m-2 s-1 KPa-1).
    psoil50 : float, optional
        Soil water potential value at 50% Lp. The default water potential is -733.6.
    Lpmax : float, optional
        Max Lp.if False, Lp = 4*Lpmin (From Henzler et al,1999 data). 
        Otherwise a value needs to be provided
        The default is False.

    Returns
    -------
    Lp : float
        Root radial conductivity (Kg m-2 s-1 KPa-1).
        
    Reference:
    ----------
    Caldeira, C. F., et al. (2014). "Circadian rhythms of hydraulic conductance 
    and growth are enhanced by drought and improve plant performance.
    " Nature Communications 5.
    
    Henzler, T., et al. (1999). "Diurnal variations in hydraulic 
    conductivity and root pressure can be correlated with the 
    expression of putative aquaporins in the roots of Lotus japonicus." 
    Planta 210(1): 50-60.

    """
    
    if not Lpmax:
        Lpmax=4*Lpmin
        
    Lp=Lpmax/(1+np.exp(-m*(np.abs(psoil)-np.abs(psoil50))))+Lpmin
    
    if Lp<Lpmin:
        Lp=Lpmin
    elif Lp>Lpmax:
        Lp=Lpmax
    
    return Lp

def varLp_flora(flux,lv,layer_depth,aroot,manage,Lpmax=0.046152,Lpmin=1.52532e-3):
    """
    Function to vary Lp based on flora pulse experiment data fitting. Ref(XX) 

    Parameters
    ----------
    flux : float. vector
        Root water uptake (kg/m2soil/s).
    lv : float
        root lenght density (m_root m-3 soil)
    layer_depth: float
        soil layer depth (m)
    aroot: float
        root radious (m)
    Lpmax : float, optional
        Maximum root hydraulic conducttivity (kg/m2 root/hour/bar). The default is 0.046152.
    Lpmin : float, optional
        Minimum root hydraulic conducttivity (kg/m2 root/hour/bar). The default is 1.52532e-3.

    Returns
    -------
    Lp : float. vector
        Actual root hydraulic conductivity (kg m-2 root s-1 kPa-1).

    """
    
    #Prevent from negative flux values
    flux=np.abs(flux)
    #Create container for the data
    Lp=np.zeros(flux.shape)
    for col in range(2): 
        for c,f in enumerate(flux[:,col]):
            
             if f==0.0: #Prevent division by zero
                flux_rb=0.0
             elif col==0: 
                 flux_rb=f/(lv[c,col]*layer_depth[c]*2*np.pi*aroot)*3600#root water flux on a root a based
                                                #(kg/m2root/hour)
             elif col ==1:
                     flux_rb=f/(lv[c,col]*layer_depth[c]*2*np.pi*aroot)*3600#root water flux on a root a based
                                                #(kg/m2root/hour)
             if flux_rb<=0.5:  
                Lp[c,col]=flux_rb/10
             else:
                Lp[c,col]=Lpmax
                
             if Lp[c,col]<Lpmin:
                Lp[c,col]=Lpmin
        
    return Lp/3600/100# root hydraulic conductivity (kg m-2_root s-1 kPa-1)
    
    


# if __name__ == '__main__':
    
#     flux=np.array([1,3,0])
#     lv=np.array([3,4,3])
#     Lp=varLp_flora(flux,lv,0.1,0.0003)
#         Tl=Tleaf(380.86,1.304,21.12,55.32,97,0.03,0.1466*1.6,gs_ad=0,epsilon_s=0.98)
         
#     import pandas as pd
#     dirfarq='C:\Trabajo\Australia\Papers\HeatWave\PythonAnalysis\Exp-Site-Soil-Comparison\exp-noirri-Barossa\input\tree'
#     ta=25
#     rint=400
#     ca=400
#     zalfa1=-33
#     zalfa2=12259821
#     LAI=1
#     Pa=101.3
#     vpd=0.5
#     potr=-1000
#     sf=0.003
#     farq=pd.read_csv(dirfarq+'\farquhar.csv',header=None,index_col=0)
#     gs,asim,psileaf,evap,ci=tuzet_joint(ta,rint,ca,zalfa1,
#                                         zalfa2,LAI,Pa,vpd,
#                                         farq.loc['g0p'],
#                                         farq.loc['mp'],
#                                         sf,potr,
#                                         farq.loc['THETA'],
#                                         farq.loc['ALFA'],
#                                         farq.loc['CTES'],
#                                         farq.loc['DELHTES'],
#                                         farq.loc['ckc'],
#                                         farq.loc['DELHKC'],
#                                         farq.loc['cKO'],
#                                         farq.loc['DELHKO'],
#                                         farq.loc['CRD'],
#                                         farq.loc['DELHRD'],
#                                         farq.loc['CVCX'],
#                                         farq.loc['DELHVCX'],
#                                         farq.loc['CJX'],
#                                         farq.loc['DELHJX'],
#                                         farq.loc['CTPU'],
#                                         farq.loc['DELHTPU'])