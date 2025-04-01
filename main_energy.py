# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:05:24 2019

@author: Omar
"""

"""
SPAC-TREE

SPAC model for tree crops.



"""
import pandas as pd
import numpy as np
import warnings
import os
import csv
import datetime
import radiation_partitioning as rad
from sun_riseSet import sun_riseSet as srs
from find_nearest import find_nearest as fn
import esoilFunc
from EToFAO import EToFAO
import tcanopy 
import modelfunc




############################ INITIALIZATION ################################
def main():
    print('Running program')
    ini_time=datetime.datetime.now()
    #Read folder name for simulations
    file=open('experiment_folder_name.txt','r')
    exp_folder=str(file.read())
    file.close()
    #dir_simulation=os.getcwd()+'/'+exp_folder+'/'
    dir_simulation=os.path.dirname(exp_folder)
    dir_spac_input=dir_simulation+'/input/'
    #import localization data
    read=csv.reader(open(dir_spac_input+'meteo/'+'coordinates.csv'))
    coord={}
    for r in read:
        coord[r[0]]=float(r[1])#coordinates
        
    
    #TREE PARAM
        
    read=csv.reader(open(dir_spac_input+'tree/'+'farquhar.csv'))
    farq={}
    for r in read:
        farq[r[0]]=float(r[1])#Farquhar param 
        
     
    read=csv.reader(open(dir_spac_input+'tree/'+'spac.csv'))
    spac={}
    for r in read:
        spac[r[0]]=float(r[1])#spac param (resistances)
    
    read=csv.reader(open(dir_spac_input+'tree/'+'treecharact.csv'))
    tree={}
    for r in read:
        tree[r[0]]=float(r[1])#tree and orchard characteristics
    
    #import irrigation options
    read=csv.reader(open(dir_spac_input+'irrigation/'+'manage.csv'))
    manage={}
    for r in read:
        manage[r[0]]=float(r[1])
    read=csv.reader(open(dir_spac_input+'irrigation/'+'irrioptions.csv'))
    irriopt={}
    for r in read:
        irriopt[r[0]]=r[1]
    
    #import soil evaporation, characteristics and initial swc
    read=csv.reader(open(dir_spac_input+'soil/'+'soil_evap_param.csv'))
    soil_evap_param={}
    for r in read:
        soil_evap_param[r[0]]=float(r[1])
    
    initswc=pd.read_csv(dir_spac_input+'soil/'+'initialswc.csv',index_col=0)#inital swc and root parameters
    
    soil=pd.read_csv(dir_spac_input+'soil/'+'soilcharact.csv')#soil layer depth and Campbell characteristics (Campbell & Norman, 1998)
    
    read=csv.reader(open(dir_spac_input+'soil/'+'soilthermal.csv'))
    soilthermal={}
    for r in read:
        soilthermal[r[0]]=float(r[1])
    
    #Soil and vegetation spectral characteristcs
    read=csv.reader(open(dir_spac_input+'radiation/'+'vegetationspectra.csv'))    
    spectraVeg={}
    for r in read:
        spectraVeg[r[0]]=float(r[1])
    
    read=csv.reader(open(dir_spac_input+'radiation/'+'soilspectra.csv'))    
    spectraGrd={}
    for r in read:
        spectraGrd[r[0]]=float(r[1])
    
    #G dsitribution in the tree
    read=csv.reader(open(dir_spac_input+'radiation/'+'inputG.csv'))
    inputG=[]
    for r in read:
        inputG.append(float(r[0]))
    inputG=np.array(inputG)
    
    #Radiation interception option (Hedgerow or single tree)    
    read=csv.reader(open(dir_spac_input+'radiation/'+'radiationInterceptionOption.csv'))
    hedgerow=[]
    for r in read:
        hedgerow.append(float(r[1]))
    hedgerow=np.array(hedgerow)
    #Weather file
    meteo=pd.read_csv(dir_spac_input+'meteo/'+'weather.csv')
    meteo['Time']=pd.to_datetime(meteo['Time'])
    
    #Check for missing data into the weather file
    meteo_check=meteo.applymap(np.isreal).all()
    if not all (meteo_check):
        pos=np.where(meteo_check==False)
        raise ValueError('The weather file contains missing values in column(s) '+ 
                         str(meteo_check.index[pos[0]]))
    
    #Budbreak dates
    budbreak=pd.read_csv(dir_spac_input+'tree/'+'budbreakdates.csv')
    budbreak['Date']=pd.to_datetime(budbreak['Date'])
    
    #Leaf area and root length density evolution
    LAImeas=pd.read_csv(dir_spac_input+'tree/'+'LAI.csv')#leaf area index
    LAImeas['Date']=pd.to_datetime(LAImeas['Date'])
    lv=pd.read_csv(dir_spac_input+'tree/'+'rld.csv')#root lenght density
    lv['Time']=pd.to_datetime(lv['Time'])
    
    #Initial root length density (Lv)
    lv_day=lv[lv['Time']==meteo.loc[0,'Time']]   
    lv_day=np.array(lv_day[['Lv_1','Lv_2']])
    
    if lv_day.size==0:
        
        raise ValueError('Initial root lenght density must match initial simulation date. '\
                         'Use the initial date in the weather file')
    ##########################LAI evolution########################################
    #create the necessary time variables
    meteo['DOY']=meteo['Time'].dt.dayofyear
    meteo['hour']=meteo['Time'].dt.hour
    meteo['minute']=meteo['Time'].dt.minute
    meteo['hour_minute']=meteo['hour']+meteo['minute']/60
    timeRecord=meteo.loc[meteo.index[1],'hour_minute']-meteo.loc[meteo.index[0],'hour_minute']#Wheather time interval
    
    #Compute the total number of days
    
    ndays=meteo.shape[0]/(1/timeRecord*24)
    
    if not ndays.is_integer():
        
        raise Exception ('The total number of days is not an integer.\n' \
                         'Check the weather file for inconsistencies')
    else:
        ndays=int(ndays)
        
       
    
    if ndays>1:
        if LAImeas.shape[0]==1:
            
            raise ValueError ('When simulating more than one day at least two LAI '\
                               'dates must be included')
            
    if LAImeas.shape[0]>1:
        
       
        
        #Create an array of LAI based on the measurements
        timeRecord=meteo.loc[meteo.index[1],'hour_minute']-meteo.loc[meteo.index[0],'hour_minute']#Wheather time interval
        datedifIni=(LAImeas.loc[0,'Date']-budbreak.loc[0,'Date'])/datetime.timedelta(days=1)
        
            
        if datedifIni>60:
            
            warnings.warn('Bud break date is 60 days before the first LAI value.'\
                         ' Reconsider to change it')
        
        
        if datedifIni<0:
            raise Exception('Budbreak date must be earlier than the first LAI measurement')
            
        LAImeas['dateDif']=LAImeas['Date'].diff()/datetime.timedelta(days=1)#days difference between measurements
        LAImeas.loc[0,'dateDif']=(LAImeas.loc[0,'Date']-budbreak.loc[0,'Date'])/datetime.timedelta(days=1)
        LAImeas['laiDif']=LAImeas['LAI'].diff()
        LAImeas.loc[0,'laiDif']=LAImeas.loc[0,'LAI']
        LAImeas['slope']=LAImeas['laiDif']/LAImeas['dateDif']#slope to calculate LAI evolution
        
        index=0
        for dif in LAImeas['dateDif']:
            ndaysLAI=np.arange(0,dif)
            if index==0:
                lai_range=ndaysLAI*LAImeas.loc[index,'slope']      
                laievol=np.array(lai_range)
            else:
                lai_range=LAImeas.loc[index-1,'LAI']+ndaysLAI*LAImeas.loc[index,'slope']        
                laievol=np.append(laievol,lai_range)    
            index=index+1
            
        laievol=np.array([laievol]).T #daily LAI evolution
        laievol=np.append(laievol,LAImeas.iloc[-1,1])
        laievol=np.repeat(laievol,1/timeRecord*24) 
        laievol=np.array([laievol]).T #daily LAI evolution 
        
        #find coindiding dates between weather and lai
        laidates=pd.date_range(LAImeas.loc[0,'Date'],periods=laievol.shape[0],
                               freq=str(round(timeRecord*60))+'min')
        laievol=pd.DataFrame(laievol,index=laidates,columns=['lai'])
        common=laidates.intersection(meteo['Time'])
        
        rows=meteo[(meteo['Time']>=str(common[0]))& (meteo['Time']<=str(common[-1]))].index.values
        
        LAI=np.zeros([meteo['DOY'].shape[0],1])
        rows=np.array(rows)
        LAI[rows]=laievol.loc[common]
        LAI=np.squeeze(LAI)
    else:
        LAI=np.repeat(LAImeas.loc[0,'LAI'],24/timeRecord)
        
    ###############################FROMAT IRRIGATION###############################
    
    
    if int(irriopt['irriopt'])==0:
        irriday=pd.read_csv(dir_spac_input+'irrigation/'+'irriFile.csv')
        irriday['Time']=pd.to_datetime(irriday['Time'])
            
    elif int(irriopt['irriopt'])==1:
        ndays
        irricycle=int(irriopt['irricycle'])
        irrini=datetime.datetime.strptime(irriopt['irriDateStart'],'%m/%d/%Y')
        irriTime=[irrini+datetime.timedelta(days=x*irricycle) for x in range(ndays)]
        irriday=pd.DataFrame({'Time':irriTime,'Irrivol':np.repeat(irriopt['irrivol'],ndays)})
    
    irriHourini=int(irriopt['irriHourStart'])
    irriHours=int(irriopt['irriHours'])
    irriHourend=int(irriHourini+irriHours)     
    
    
    ###########################PARAMETER CALCULATION###############################
    
    #Residual soil water content to stop flux extraction from the roots
    
    theta_res=modelfunc.theta_residual(np.array(soil['pe']),
                                        np.array(soil['b']),
                                        np.array(soil['sat']))
    
    #Intial stem water potential for irrigation triggering
    pstem=0.0#tree['initial_pstem']
    
    #Initial swc
    swc=np.array(initswc)
    
    #soil layer depth
    layer_depth=np.array(soil['layer_depth'])
    
    #Calculate atmospheric pressure
    Pa=101.3*np.exp(-coord['altitude']/8200) 
    
    #Initial values for soil evaporation calculation
    CEwet=soil_evap_param['Cewet']#If zero, soil is wetted
    CEdry=soil_evap_param['Cedry']#If zero, soil is wetted
    
    twet=0
    tdry=0
    
    #fraction of ground cover
    f_c=tree['canopywidth']/tree['rowdist']
    
    # Canopy width to height ratio
    w_c_ratio = tree['canopywidth'] / (tree['htop'] - tree['hbase']) 
    
    #Coordinates
    lat= coord['lat']#Latitude
    lon=coord['long']#longitude
    stdlon=coord['standardLong']
    londiff=lon-stdlon#standard meridian difference
    
    #Anemomter height
    aneHeight=coord['aneHeight']
    
    #vegetation indeces
    x_LAD=tree['x_LAD']#Canopy angle distribution 
    
    ##create indexes for slicing day by day
    drange=np.arange(0,meteo.shape[0],int(24/timeRecord))
    if not [drange]:
        drange=[0]
    day=0 #Initial day 
    
    ###############################CONTAINER FOR THE DATA##########################
    
    
    Rn_sw_veg_soil=np.zeros([meteo['DOY'].shape[0],2]) #vegetation and soil shortwave Rn
    Rn_par_dir_dif=np.zeros([meteo['DOY'].shape[0],2]) #vegetation PAR radiation
    Tc=np.zeros(meteo['DOY'].shape[0])#canopy temperature
    Tleaf=np.zeros([meteo['DOY'].shape[0],2])
    Ep=np.zeros(meteo['DOY'].shape[0])#Plant transpiration soil surface based
    Ep_plant=np.zeros(meteo['DOY'].shape[0])#Instantaneous plant transpiration l/h
    VPD=np.zeros(meteo['DOY'].shape[0])
    Pstem=np.zeros(meteo['DOY'].shape[0])#collar water potential
    Pleaf_sun=np.zeros(meteo['DOY'].shape[0])#sunlit leaf water potential
    Pleaf_shade=np.zeros(meteo['DOY'].shape[0])#shade leaf water potential
    Gc=np.zeros(meteo['DOY'].shape[0])#sunlit leaves conductance
    Gs_sun=np.zeros(meteo['DOY'].shape[0])#shaded leaves conductance
    Gs_shade=np.zeros(meteo['DOY'].shape[0])#canopy conductance
    Pn=np.zeros(meteo['DOY'].shape[0])#canopy photosynthesis
    SWCwet=np.zeros([meteo['DOY'].shape[0],swc.shape[0]])
    SWCdry=np.zeros([meteo['DOY'].shape[0],swc.shape[0]])
    ksoil=np.zeros(swc.shape)
    psoil=np.zeros(swc.shape)
    rsoil=np.zeros(swc.shape)
    rr=np.zeros(swc.shape)
    flux=np.zeros(swc.shape)
    Lp=np.zeros(swc.shape)
    LpMean=np.zeros(meteo['DOY'].shape[0])
    flux_tot=np.zeros(meteo['DOY'].shape[0])
    
    #Daily
    Eswet=np.zeros(ndays)
    Esdry=np.zeros(ndays)
    ETo=np.zeros(ndays)
    Qd=np.zeros(ndays)
    IrriVol=np.zeros(ndays)
    RainDay=np.zeros(ndays)
    #######################DAILY CALCULATIONS######################################
    
    for d in drange: 
        dend=int(d+24/timeRecord)
        meteo_day=meteo.iloc[d:dend,:]
        
        ################## SOLAR RADIATON ####################################
        # Compute illumination angles
        sza, saa = rad.calc_sun_angles(np.full(meteo_day['DOY'].shape,
                                           lat),
                                   np.full(meteo_day['DOY'].shape,
                                           lon),
                                   np.full(meteo_day['DOY'].shape,
                                           stdlon),
                                   meteo_day['DOY'],
                                   meteo_day['hour_minute'])
                                   # Relative angle row-wind
        
        
        Sdn =np.array( meteo_day['Srad'] ) # Solar irradiance (W m-2)
        p = 10 * np.full(meteo_day['DOY'].shape,Pa)  # Atmospheric pressure in hPa
        
        # Avoid unrealistic values at sundown and sunset
        Sdn[sza >= 90] = 0
        
        # Estimates the direct and diffuse solar radiation
        difvis, difnir, fvis, fnir = rad.calc_difuse_ratio(Sdn, sza, press=p)
        # Solar diffuse irradiance ratio
        skyl = fvis * difvis + fnir * difnir
        Sdn_dir = (1. - skyl) * Sdn
        Sdn_dif = skyl * Sdn
        
        #Estimates direc and diffuse in the PAR waveband
        skylPAR = fvis * difvis
        Sdn_dirPAR = (1. - skylPAR) * Sdn
        Sdn_difPAR = skylPAR * Sdn
        
        
        ###################  RADIATION INTERCEPTION ###########################
        
        if LAI[d]>0:
            
            
            #Calculate sunrise and sunset and find rows
            sunrise,sunset=srs(meteo_day.loc[meteo_day.first_valid_index(),'DOY'],
                               lat,
                               londiff)
            
            sunrise=fn(meteo_day['hour_minute'], sunrise)
            sunset=fn(meteo_day['hour_minute'], sunset)
            index_sunrise=np.where(meteo_day['hour_minute']==sunrise)[0][0]
            index_sunset=np.where(meteo_day['hour_minute']==sunset)[0][0]
            rows=meteo_day.index.values
            
            # Compute illumination angles
            # sza, saa = rad.calc_sun_angles(np.full(meteo_day['DOY'].shape,
            #                                    lat),
            #                            np.full(meteo_day['DOY'].shape,
            #                                    lon),
            #                            np.full(meteo_day['DOY'].shape,
            #                                    stdlon),
            #                            meteo_day['DOY'],
            #                            meteo_day['hour_minute'])
            #                            # Relative angle row-wind
            if hedgerow:
    #####################RADIATION HEDGEROW########################################            
                # Use row clumping index
                psi = tree['rowdirec'] - saa
                
                # Row crops model
                Omega = rad.calc_omega_rows(LAI[rows],
                                                f_c,
                                                theta=sza,
                                                psi=psi,
                                                D=w_c_ratio,
                                                x_LAD=x_LAD,
                                                isLAIeff=False)
                
               
                
                
                Sdn =np.array( meteo_day['Srad'] ) # Solar irradiance (W m-2)
                p = 10 * np.full(meteo_day['DOY'].shape,Pa)  # Atmospheric pressure in hPa
                
                # Avoid unrealistic values at sundown
                Sdn[sza >= 90] = 0
                Omega[sza > 90] = 1
                sza[sza > 90] = 90
                Omega[~np.isfinite(LAI[rows])] = np.nan
                Omega[np.array(meteo_day['Srad']) < 0] = 1
                
                #Avoid calculation errors at low radiation
                if np.any(Omega==0):
                
                    fake_Omega_pos=np.where(Omega==0.0)
                    Omega[fake_Omega_pos[0]]=1
            
                #Local LAI 
                F=LAI[rows]/f_c
                
                # Calculate the effective LAI at illumination angle
                
                if tree['LAI_eff']:
                    
                    LAIeff = F * Omega
                else:
                    LAIeff=None
                # # Estimates the direct and diffuse solar radiation
                # difvis, difnir, fvis, fnir = rad.calc_difuse_ratio(Sdn, sza, press=p)
                # # Solar diffuse irradiance ratio
                # skyl = fvis * difvis + fnir * difnir
                # Sdn_dir = (1. - skyl) * Sdn
                # Sdn_dif = skyl * Sdn
                
                # #Estimates direc and diffuse in the PAR waveband
                # skylPAR = fvis * difvis
                # Sdn_dirPAR = (1. - skylPAR) * Sdn
                # Sdn_difPAR = skylPAR * Sdn
                
                
                #Net short wave radiation 
                Rn_sw_veg_PT,Rn_dif_veg_PT, Rn_sw_soil_PT, Rn_par_dir,Rn_par_dif = rad.calc_Sn_Campbell(
                                                    LAI[rows],
                                                    sza,
                                                    Sdn_dir,
                                                    Sdn_dif,
                                                    fvis,
                                                    fnir,
                                                    np.full(LAI[rows].shape,
                                                            spectraVeg['rho_leaf_vis']),
                                                    np.full(LAI[rows].shape,
                                                            spectraVeg['tau_leaf_vis']),
                                                    np.full(LAI[rows].shape,
                                                            spectraVeg['rho_leaf_nir']),
                                                    np.full(LAI[rows].shape,
                                                            spectraVeg['tau_leaf_nir']),
                                                    np.full(LAI[rows].shape,
                                                            spectraGrd['rsoilv']),
                                                    np.full(LAI[rows].shape,
                                                            spectraGrd['rsoiln']),
                                                    x_LAD=x_LAD,
                                                    LAI_eff=LAIeff
                                                        )
                Rn_sw_veg_PT[0:index_sunrise]=0
                Rn_sw_veg_PT[index_sunset:]=0
                
                Rn_dif_veg_PT[0:index_sunrise]=0
                Rn_dif_veg_PT[index_sunset:]=0
                
                Rn_sw_soil_PT[0:index_sunrise]=0
                Rn_sw_soil_PT[index_sunset:]=0
                Rn_sw_veg_soil[rows,:]=np.array([Rn_sw_veg_PT,Rn_sw_soil_PT]).T# W m-2 soil
                
                Rn_par_dir[0:index_sunrise]=0
                Rn_par_dir[index_sunset:]=0
                Rn_par_dif[0:index_sunrise]=0
                Rn_par_dif[index_sunset:]=0
                Rn_par_dir_dif[rows,:]=np.array([Rn_par_dir,Rn_par_dif]).T# W PAR m-2 soil
                
            else:
    ###############RADIATION SPHERICAL TREE CANOPY#################################
                
                Sdn =np.array( meteo_day['Srad'] ) # Solar irradiance (W m-2)
                p = 10 * np.full(meteo_day['DOY'].shape,Pa)  # Atmospheric pressure in hPa
                
                # Avoid unrealistic values at sundown and sunset
                Sdn[sza >= 90] = 0
                #Container for daily values
                Rn_sw_veg_PT=np.zeros(Sdn.shape[0])
                Rn_dif_veg_PT=np.zeros(Sdn.shape[0])
                Rn_sw_soil_PT=np.zeros(Sdn.shape[0])
                Rn_par_dir=np.zeros(Sdn.shape[0])
                Rn_par_dif=np.zeros(Sdn.shape[0])
                LAIsun=np.zeros(Sdn.shape[0])
                LAIshade=np.zeros(Sdn.shape[0])
                
                # # Estimates the direct and diffuse solar radiation
                # difvis, difnir, fvis, fnir = rad.calc_difuse_ratio(Sdn, sza, press=p)
                # # Solar diffuse irradiance ratio PAR+NIR
                # skyl = fvis * difvis + fnir * difnir
                # Sdn_dir = (1. - skyl) * Sdn
                # Sdn_dif = skyl * Sdn
                
                #Compute radiation intercepted between sunrise and sunset
                #index_rad=index_sunrise
                
                for s in sza[sza<=90]:
                    index_rad=np.where(s==sza)
                    index_rad=index_rad[0]
                    if Sdn_dirPAR[index_rad]>0:
                        rad_result_NIR=rad.rint_spheroid(s,tree['rowdist'],tree['plantspace'],
                                          tree['ratio_rx_rz'],
                                          tree['htop'],tree['hbase'],
                                          inputG,LAI[rows[0]],
                                          Sdn_dir[index_rad],Sdn_dif[index_rad])
                    else:
                        rad_result_NIR=np.zeros(2) 
                    Rn_sw_veg_PT[index_rad]=rad_result_NIR[0]+rad_result_NIR[1]# W m-2 soil
                    Rn_dif_veg_PT[index_rad]=rad_result_NIR[1]# W m-2 soil
                    Rn_sw_soil_PT[index_rad]=Sdn[index_rad]-Rn_sw_veg_PT[index_rad]# W m-2 soil
                    
                    #index_rad=index_rad+1
                    
                
                Rn_sw_veg_PT[0:index_sunrise]=0
                Rn_sw_veg_PT[index_sunset:]=0
                
                Rn_dif_veg_PT[0:index_sunrise]=0
                Rn_dif_veg_PT[index_sunset:]=0
                
                Rn_sw_soil_PT[0:index_sunrise]=0
                Rn_sw_soil_PT[index_sunset:]=0
                Rn_sw_veg_soil[rows,:]=np.array([Rn_sw_veg_PT,Rn_sw_soil_PT]).T
                
                # # Solar diffuse irradiance ratio PAR
                # skylPAR = fvis * difvis
                # Sdn_dirPAR = (1. - skylPAR) * Sdn
                # Sdn_difPAR = skylPAR * Sdn
                
                #index_rad=index_sunrise
                for s in sza[sza<=90]:
                    index_rad=np.where(s==sza)
                    index_rad=index_rad[0]
                    
                    if Sdn_dirPAR[index_rad]>0:
                        APARDIR,APARDIF,laisun,laishade=rad.rint_spheroid(s,tree['rowdist'],tree['plantspace'],
                                          tree['ratio_rx_rz'],
                                          tree['htop'],tree['hbase'],
                                          inputG,LAI[rows[0]],
                                          Sdn_dirPAR[index_rad],Sdn_difPAR[index_rad])
                    else:
                        APARDIR=0
                        APARDIF=0
                        laisun=0
                        laishade=LAI[d]
                    Rn_par_dir[index_rad]=APARDIR# W PAR m-2 soil
                    Rn_par_dif[index_rad]=APARDIF# W PAR m-2 soil
                    LAIsun[index_rad]=laisun
                    LAIshade[index_rad]=laishade
                    #index_rad=index_rad+1
                    
                Rn_par_dir[0:index_sunrise]=0
                Rn_par_dir[index_sunset:]=0
                Rn_par_dif[0:index_sunrise]=0
                Rn_par_dif[index_sunset:]=0
                
                LAIsun[0:index_sunrise]=0
                LAIsun[index_sunset:]=0
                LAIshade[0:index_sunrise]=0
                LAIshade[index_sunset:]=0
                
                Rn_par_dir_dif[rows,:]=np.array([Rn_par_dir,Rn_par_dif]).T
            
                
        else:
            
            Rn_sw_veg_PT=np.zeros(np.int(24/timeRecord))
            Rn_sw_soil_PT=np.zeros(np.int(24/timeRecord))
            Rn_sw_veg_soil[rows,:]=0.0
            
            Rn_par_dir=np.zeros(np.int(24/timeRecord))
            Rn_par_dif=np.zeros(np.int(24/timeRecord))
            Rn_par_dir_dif[rows,:]=0
            
        ####################SOIL TEMPERATURE EVOLUTION#########################
        
        #Compute soil sensible heat evolution
        t_G=np.arange(0,86400,timeRecord*3600)#time array for G calculations
        G=tcanopy.G_surface(meteo_day['Temp'],
                             soilthermal['kc'],
                             soilthermal['kd'],
                             t_G,
                             t0=8)# (W m-2)
            
        ####################IRRIGATION#########################################
        
        if int(irriopt['irriopt'])==2:
            #Trigger irrigation when the minimum stem water potential of the
            #previuos day passes the irriation threshold
            if day>0:
                d_previous=d-np.int(24/timeRecord)
                if np.min(Pstem[d_previous:d])<=np.float(irriopt['irriThres']):
                
                    irrivol=np.float(irriopt['irrivol'])
            
                else:
                
                    irrivol=0.0
                    
            else:
                
                irrivol=0.0
            
        
        else:
            
            if any(meteo_day.loc[d,'Time']==irriday['Time']):
                
                index=irriday[irriday['Time']==meteo_day.loc[d,'Time']]
                irrivol=np.float(index['Irrivol'])
            else:
                irrivol=0.0
                
        IrriVol[day]=irrivol
        RainDay[day]=np.sum(meteo_day['Rain'])
        ##################ETo################################################## 
        Ta=meteo_day['Temp']# Air temperature (C)
        Hr=meteo_day['Hr']#Relative humidity (%)
        u=meteo_day['Wind']#wind velocity (m s-1)
        u_htop=tcanopy.u_topcanopy(u,aneHeight,tree['htop'])
        #Force wind at the top of teh canopy to avoid unrealistic aerodynamic resistances
        if any(u_htop<=0.8):
            u_htop[u_htop<=0.8]=0.8
            
        Rs=meteo_day['Srad']#instantaneous solar radiation (W m-2)
        ETo[day]=EToFAO(Ta,u,Hr,Rs,meteo_day.loc[d,'DOY'],lat,timeRecord,AneHeight=aneHeight)     
        #################SOIL EVAPORATION######################################
        
            
         
        vpd=esoilFunc.es(Ta)*(1-Hr/100)#vapor pressure deficit (kPa)
        VPD[d:dend]=vpd
        Rn_soil_day=(np.sum(Rn_sw_soil_PT)*3600*timeRecord)/1e6 #Daily net radiation 
                                                                #reaching the soil
        #Rn_veg_day=(np.sum(Rn_sw_veg_PT)*3600*timeRecord)/1e6 #Daily net radiation 
                                                                #reaching the canopy
        #Rs_day =(np.sum(Rs)*3600*timeRecord)/1e6# daily solar radiation (MJ m-2 day-1)
        Rn_veg_day_PAR=(np.sum(Rn_par_dir+Rn_par_dif)*3600*timeRecord)/1e6
        
        Rs_day_PAR=(np.sum(Sdn_dirPAR+Sdn_difPAR)*3600*timeRecord)/1e6
        
        Qd[day]=Rn_veg_day_PAR/Rs_day_PAR#Average fraction of intercepted radiation
        #Rn_soil_day=(np.sum(Rn_sw_soil_PT)*3600*timeRecord)/1e6 #Daily net radiation 
                                                                #reaching the soil
        Rn_soil_day=Rs_day_PAR*(1-Qd[day])
        
        if np.sum(meteo_day['Rain'])>=soil_evap_param['U']:
            CEwet=0
            CEdry=0
            twet=0
            tdry=0
        
        
        if irrivol>ETo[day]*(1-Qd[day])*manage['fwet']:
            CEwet=0
            twet=0
            
        
        #Initialize the time to compute second phase evaporation
        if CEwet>=soil_evap_param['U']:
            twet=twet+1
        
        if CEdry>=soil_evap_param['U']:
            tdry=tdry+1
        Eswet[day],Esdry[day],CEwet,CEdry=esoilFunc.Esoil(ETo[day],np.mean(Ta),
                                                         Pa,Qd[day],
                                                         Rn_soil_day,
                                                         np.mean(vpd),
                                                         np.mean(u),
                                                         manage['fwet'],
                                                         manage['fdry'],
                                                         CEwet,CEdry,
                                                         twet,tdry,
                                                         Ue=soil_evap_param['U'],
                                                         ce=soil_evap_param['ce'])
        
        
     
            
        #Evaporation layer to limit Es when the soil is dry. 
        
        evapLayer=swc[0,:]*layer_depth[0]*1000#(mm)
        #Limit evaporation to the available water 
        evapLayer[0]=evapLayer[0]-Eswet[day]#(mm)
        evapLayer[1]=evapLayer[1]-Esdry[day]#(mm)                                                                                          
        #Limit the minimum swc in the first layer due to evaporation and makes 
        #Esoil=0.1 residual evaporation
        if np.any(evapLayer<=0.0):
            
            dry_layer=np.where(evapLayer<=0.0)
            
            if np.all(dry_layer[0]==0):
                Eswet[day]=0.1
            elif np.all(dry_layer[0]==1):
                Esdry[day]=0.1
            else:
                Eswet[day]=0.1
                Esdry[day]=0.1
        #Correct for the Es value when changing from phase I to phase II
        if day>0:
            if twet==1:
                if Eswet[day-1]<soil_evap_param['ce']:
                    Eswet[day]=Eswet[day-1]
            if twet>1:
                if Eswet[day]>Eswet[day-1]:
                    Eswet[day]=Eswet[day-1]
            if tdry==1:
                if Esdry[day-1]<soil_evap_param['ce']:
                    Esdry[day]=Esdry[day-1]
            if tdry>1:
                if Esdry[day]>Esdry[day-1]:
                    Esdry[day]=Esdry[day-1]
                    
    ###########UPDATE ROOT LENGTH DENSITY##########################################
        
        if (lv['Time']==meteo_day.loc[d,'Time']).any(): 
            lv_day=lv[lv['Time']==meteo_day.loc[d,'Time']] 
            lv_day=np.array(lv_day[['Lv_1','Lv_2']])
        #################SUB-DAY CALCULATIONS##################################
        
        time=0
        for t in np.arange(d,dend):
            
            #Compute Rain intercepted by the canopy
            if np.sum(meteo_day['Rain'])>=1:
                if LAI[t]==0:
                    rainWet=meteo_day.loc[t,'Rain']
                    rain=meteo_day.loc[t,'Rain']
                else:
                    rainWet=modelfunc.rainfall_part(LAI[t],f_c,meteo_day.loc[t,'Rain'])
                    rain=meteo_day.loc[t,'Rain']
            else:
                rain=0.0
                rainWet=0.0
            
            #Calculate instantanoeus evaporation during the day 
            sunnytime=meteo_day[meteo_day['Srad']>10].index.values
            if t in sunnytime:
                Eswet_inst=Eswet[day]/sunnytime.shape[0]
                swc[0,0]=swc[0,0]-Eswet_inst/layer_depth[0]/1000
               
               
                Esdry_inst=Esdry[day]/sunnytime.shape[0] 
                swc[0,1]=swc[0,1]-Esdry_inst/layer_depth[0]/1000
            
            if meteo_day.loc[t,'hour']>=irriHourini and meteo_day.loc[t,'hour']<irriHourend:
                
                irrivolinstant=irrivol/(irriHours/timeRecord)#Instantaneous irrigation
                
                
                swc[0,0]=swc[0,0]+(irrivolinstant+rainWet)/layer_depth[0]/1000
                swc[0,1]=swc[0,1]+(rain)/layer_depth[0]/1000
            else:
                swc[0,0]=swc[0,0]+(rainWet)/layer_depth[0]/1000
                swc[0,1]=swc[0,1]+(rain)/layer_depth[0]/1000
                
            #Prevent for swc to be above saturation values
            if any (swc[0,:]>soil['sat'][0]):
                pos=np.where(swc[0,:]>soil['sat'][0])[0]
                swc[0,pos]=soil['sat'][0]
                
            ###Variable Lp#########
            if day == 0 and t == 0:# minimum Lp at the beginning 
                Lp[:]=spac['Lp']
                
            else:
                # Lp=modelfunc.varLp_flora(flux, 
                #                           lv_day, 
                #                           soil['layer_depth'], 
                #                           tree['aroot'],
                #                           manage) 
                # m=0.01#Slope of Lp change
                # Lpmin=4.237e-9
                # psoil50=-700
                # Lpmax=1.282e-7
                # Lp[:]=modelfunc.varLp(Pstem[t-1],m, Lpmin=Lpmin,psoil50=psoil50,Lpmax=Lpmax)
                
                
                Lp[:]=spac['Lp']
                
            LpMean[t]=np.mean(Lp)
            rrT=1/(Lp*2*np.pi*tree['aroot'])
            
            for col in range(2):
                
                psoil[:,col]=modelfunc.Psoil(swc[:,col],soil['sat'],soil['pe'],soil['b'])
                ksoil[:,col]=modelfunc.Ksoil(soil['pe'],soil['b'],psoil[:,col],soil['ks'])
                
                if any (swc[:,col]>soil['sat']):
                    pos=np.where(swc[:,col]>soil['sat'])[0]
                    psoil[pos,col]=soil.loc[pos,'pe']
                    ksoil[pos,col]=soil.loc[pos,'ks']
                
                #Compute soil resistance (with roots) for each soil layer (s m2 kPa kg-1)
                rsoil[:,col]=modelfunc.Rsoil(lv_day[:,col],tree['aroot'],\
                             ksoil[:,col],soil['layer_depth'])
                
                #Compute minimum root resistance (s m kPa kg-1) 
                # if spac['Lp']>0:
                #     m=0.01#Slope of Lp change
                #     Lpmin=4.237e-9
                #     psoil50=-700
                #     Lpmax=1.282e-7
                #     Lp=modelfunc.varLp(Pstem[t-1],m, Lpmin=Lpmin,psoil50=psoil50,Lpmax=Lpmax)
                    
                  
                #     rrT=1/(Lp*2*np.pi*tree['aroot'])
                # else:
                    
                #     rrT=modelfunc.RrootT(Ta[t],tree['aroot'],
                #                       spac['y0'],spac['a'],spac['b'])
                
                #Compute root resistance as a function of soil water content (s m kPa kg-1)
                rr[:,col]=modelfunc.RrootTheta(rrT[:,col],swc[:,col],soil['sat'],
                                        spac['alfa'],spac['beta'],spac['sigma'])
                
            #Weighted soil and root resistance in each soil compartment (s m2 kPa kg-1)   
            rsoil[:,0]=rsoil[:,0]/manage['fwet']
            rsoil[:,1]=rsoil[:,1]/manage['fdry']
            rr[:,0]=rr[:,0]/(lv_day[:,0]*manage['fwet']*soil['layer_depth'])
            rr[:,1]=rr[:,1]/(lv_day[:,1]*manage['fdry']*soil['layer_depth'])
            
            
            s1=np.sum(np.sum(psoil/(rsoil+rr)))
            s2=np.sum(np.sum(1/(rsoil+rr)))
            #xylem resistance as a function of past stem water potential(s m2 kPa kg-1)
            
            
            pstem_old=Pstem[t-1]
            
            
                
                
            k_xyl=modelfunc.k_xylem(pstem_old,Kmax=spac['Kmax'],B=spac['B'],C=spac['C'])
            k_xyl=k_xyl*tree['sapwoodarea']/tree['shootheight']
           
                
            rrl=1/k_xyl
            
            
            
            zalfa1=s1/s2
            zalfa2=(rrl+1/s2)
            
            #Weight canopy resistances as a function of sunlit and shade fractions
            if LAI[d]>0:
                
                
                    
                #Compute stomtal conductance, transpiration, water potential and CO2 at
                #leaf level
                
                if Rn_par_dir[time]>20:
                    if hedgerow:
                        
                        LAI_sun,LAI_shade=rad.LAI_sun_shade(LAI[d],Omega[time],theta=sza[time],x_LAD=1)
                       
                        par_sun_leaf=Rn_par_dir[time] + Rn_par_dif[time]#W PAR m-2 leaf
                        par_shade_leaf=Rn_par_dif[time]#W PAR m-2 leaf
                    else:
                        #Take LAI sun for the corresponding hour
                        LAI_sun=LAIsun[time]
                        LAI_shade=LAIshade[time]
                     
                        par_sun_leaf=Rn_par_dir[time] + Rn_par_dif[time]#W PAR m-2 leaf
                        par_shade_leaf=Rn_par_dif[time]#W PAR m-2 leaf
                        
                    fraction_sun=LAI_sun/LAI[d]
                
                    if fraction_sun>0:
                        
                        zalfa2_sun=zalfa2/fraction_sun
                        zalfa2_shade=zalfa2/(1-fraction_sun)
                        
                        rrl_sun=rrl/fraction_sun
                        rrl_shade=rrl/(1-fraction_sun)
                        
                    else:
                        
                        zalfa2_sun=zalfa2
                        zalfa2_shade=zalfa2
                        
                        rrl_sun=rrl
                        rrl_shade=rrl
                        
                    #Calculations from sunlit leaves
                    gs_sun,asim_sun,psileaf_sun,\
                    evap_sun,ci_sun,tleaf_sun= modelfunc.tuzet_energy(Ta[t],
                                                        par_sun_leaf,
                                                        meteo_day.loc[t,'CO2'],
                                                        u_htop[t],
                                                        Hr[t],
                                                        coord['altitude'],
                                                        tree['leaf_width'],
                                                        zalfa1,
                                                        zalfa2_sun,
                                                        LAI_sun,
                                                        farq['g0p'],
                                                        farq['mp'],
                                                        spac['sf'],
                                                        spac['potr'],
                                                        farq['THETA'],
                                                        farq['ALFA'],
                                                        farq['CTES'],
                                                        farq['DELHTES'],
                                                        farq['ckc'],
                                                        farq['DELHKC'],
                                                        farq['cKO'],
                                                        farq['DELHKO'],
                                                        farq['CRD'],
                                                        farq['DELHRD'],
                                                        farq['CVCX'],
                                                        farq['DELHVCX'],
                                                        farq['CJX'],
                                                        farq['DELHJX'],
                                                        farq['CTPU'],
                                                        farq['DELHTPU'])
                   
                    Tleaf[t,0]=tleaf_sun        
                    
                    #Calculations for shaded leaves
                    
                    gs_shade,asim_shade,psileaf_shade,\
                    evap_shade,ci_shade,tleaf_shade= modelfunc.tuzet_energy(Ta[t],
                                            par_shade_leaf,
                                            meteo_day.loc[t,'CO2'],
                                            u_htop[t],
                                            Hr[t],
                                            coord['altitude'],
                                            tree['leaf_width'],
                                            zalfa1,
                                            zalfa2_shade,
                                            LAI_shade,
                                            farq['g0p'],
                                            farq['mp'],
                                            spac['sf'],
                                            spac['potr'],
                                            farq['THETA'],
                                            farq['ALFA'],
                                            farq['CTES'],
                                            farq['DELHTES'],
                                            farq['ckc'],
                                            farq['DELHKC'],
                                            farq['cKO'],
                                            farq['DELHKO'],
                                            farq['CRD'],
                                            farq['DELHRD'],
                                            farq['CVCX'],
                                            farq['DELHVCX'],
                                            farq['CJX'],
                                            farq['DELHJX'],
                                            farq['CTPU'],
                                            farq['DELHTPU'])
                    
                    Tleaf[t,1]=tleaf_shade 
                    
                
                else:
                    
                    gs_sun=0.0
                    gs_shade=farq['g0p']
                    
                    evap_sun=0.0
                    evap_shade=1.6*farq['g0p']*vpd[t]/103*18 #transpiration(g/m2/s)
                    
                    rrl_sun=0.0
                    #rrl_shade=rrl
                    
                    LAI_sun=0.0
                    fraction_sun=0.0
                    LAI_shade=LAI[d]
                    
                    
                    asim_sun=0.0
                    asim_shade=0.0
                    
                    
                    
                    Tleaf[t,:]=modelfunc.Tleaf(Rn_sw_veg_PT[time],
                                            u_htop[t],
                                            Ta[t],
                                            Hr[t],
                                            aneHeight,
                                            tree['leaf_width'],
                                            gs_shade)
                        
                #Calculate average canopy conductance and canopy temperature
                Gs_sun[t]=gs_sun
                Gs_shade[t]=gs_shade
                Gc[t]=gs_sun*fraction_sun+gs_shade*(1-fraction_sun)
                Gc_tc=Gc[t]*LAI[d]#Gc in m2 of soil
                Tc[t]=Tleaf[t,0]*fraction_sun+Tleaf[t,1]*(1-fraction_sun)
                
                #     Tc[t]=tcanopy.tcanopy(Rn_sw_veg_PT[time],
                #                           Ta[t],
                #                           u_htop[t],
                #                           Hr[t],
                #                           Pa,
                #                           tree['htop'],
                #                           Gc_tc,
                #                           G[time])
                
                    
                
                #Calculate whole canopy assimilation
                Pn[t]=(asim_sun*LAI_sun+asim_shade*LAI_shade)*tree['plantspace']*tree['rowdist']
          
            
                    
            
                #Plant transpiration (kg/m2 soil/s)
                Ep[t]=(evap_sun*LAI_sun+evap_shade*LAI_shade)/1000
                
                #Instantaneous Plant Transpiration l/h
                Ep_plant[t]=(evap_sun*LAI_sun*tree['plantspace']*tree['rowdist'] 
                            +evap_shade*LAI_shade*tree['plantspace']*tree['rowdist'])*3.6
                
                #Collar water potential
                
                pcollar=(s1-Ep[t])/s2
                #Stem water potential
                pstem=pcollar-Ep[t]*rrl
                Pstem[t]=pstem
                
                #compute sun and shaded leaf water potential
                # pleaf_shade=pstem-(evap_shade*LAI_shade/1000)*rrl_shade#Shaded Leaf potential (kPa)
                if rrl_sun==0.0:
                     psileaf_shade=pstem
                     psileaf_sun=psileaf_shade
                # else:             
                #     pleaf_sun=pstem-(evap_sun*LAI_sun/1000)*rrl_sun#Sunlit Leaf potential (kPa)
                
                  
                Pleaf_sun[t]=psileaf_sun
                Pleaf_shade[t]=psileaf_shade
                
                
                # Calculate water extracted from the root
                
                for col in range(2):
                    
                    flux[:,col]=(psoil[:,col]-pcollar)/(rsoil[:,col]+rr[:,col])#Root water extraction (kg s-1 m-2)
                    
                    swc[:,col]=swc[:,col]-flux[:,col]*timeRecord*3600/1000/soil.loc[:,'layer_depth']
                    #Equal flux to zero when the soil water content is below its residual value
                    if np.any(np.less_equal(swc[:,col],theta_res)):
                        
                        row_theta_res=np.where(np.less_equal(swc[:,col],theta_res))
                        swc[row_theta_res,col]=theta_res[row_theta_res]
                    
                flux_tot[t]=np.sum(flux[:,0])+np.sum(flux[:,1])
      #############################CALCULATE WATER MOVEMENT THROUGH LAYERS#########              
            
            #Update of soil water potential & conductivity between layers    
            for col in range(2):
                    
                    for layer in range(swc.shape[0]-1):
    
                        psoilb = modelfunc.Psoil(swc[layer+1,col],
                                                 soil.loc[layer+1,'sat'],
                                                 soil.loc[layer+1,'pe'],
                                                 soil.loc[layer+1,'b'])#matric potential lower layer
            
                        kb = modelfunc.Ksoil(soil.loc[layer+1,'pe'],
                                                   soil.loc[layer+1,'b'],
                                                   psoilb,
                                                   soil.loc[layer+1,'ks'])#conductivity lower layer
                
            
                        psoila = modelfunc.Psoil(swc[layer,col],
                                                 soil.loc[layer,'sat'],
                                                 soil.loc[layer,'pe'],
                                                 soil.loc[layer,'b'])#matric potential upper
            
                        ka = modelfunc.Ksoil(soil.loc[layer,'pe'],
                                                   soil.loc[layer,'b'],
                                                   psoila,
                                                   soil.loc[layer,'ks'])#conductivity upper layer
                        
                        
                        if kb > soil.loc[layer+1,'ks']:
                            kb = soil.loc[layer+1,'ks']
                            psoilb=soil.loc[layer+1,'pe']
                        
            
                        if ka > soil.loc[layer,'ks']:
                            ka = soil.loc[layer,'ks']
                            psoila=soil.loc[layer,'pe']
                        
            
                        
                        
            
                        kmed = 0.5 * (ka + kb)
            
                        fluxun=kmed*(psoila-psoilb)+9.81*kmed #Flux between layers (kg s-1 m-2)
                        
                        #Amount of water transferred at the corresponding period  (m3/m3)
                        #water_tansf=fluxun*3600*timeRecord/layer_depth[layer]/1000 
                        water_tansf=fluxun*3600*timeRecord/1000 
                        swc[layer+1, col] = swc[layer+1, col] + water_tansf/layer_depth[layer+1]
                        swc[layer, col] = swc[layer, col] - water_tansf/layer_depth[layer]
                         
                               
                    if swc[-1,col]>soil.loc[soil.shape[0]-1,'fc']:
                         dp=soil.loc[soil.shape[0]-1,'swcon']*(swc[soil.shape[0]-1,col]-soil.loc[soil.shape[0]-1,'fc'])#deep percolation
                         swc[-1,col]=swc[-1,col]-dp
                                
                            
                     
            SWCwet[t]=swc.T[0]
            SWCdry[t]=swc.T[1]
            #Update soil water potential and conductivity
            for col in range(2):
                
                  psoil[:,col]=modelfunc.Psoil(swc[:,col],
                                               soil['sat'],
                                               soil['pe'],
                                               soil['b'])
                  
                  ksoil[:,col]=modelfunc.Ksoil(soil['pe'],
                                               soil['b'],
                                               psoil[:,col],
                                               soil['ks']) 
                  
            
            
            time=time+1
            
            
        day=day+1
        print('day: '+str(day))
    
    #Format data to save it
    output={'Time':meteo['Time'],
            'Rint_veg (W m-2)':Rn_sw_veg_soil[:,0], 
            'Rint_soil (W m-2)':Rn_sw_veg_soil[:,1], 
            'Rint_dir(W m-2)':Rn_par_dir_dif[:,0],
            'Rint_dif (W m-2)':Rn_par_dir_dif[:,1],
            'Tleaf_sunlit(C)':Tleaf[:,0],
            'Tleaf_shade(C)':Tleaf[:,1],
            'Tcanopy':Tc,
            'Ep (l h-1)':Ep_plant,
            'VPD (kPa)':VPD,
            'Pstem (kPa)':Pstem,
            'Pleaf_sun (kPa)':Pleaf_sun,
            'Pleaf_shade (kPa)':Pleaf_shade,
            'Gc (mol m-2 s-1)':Gc,
            'Gs_sun (mol m-2 s-1)':Gs_sun,
            'Gs_shade (mol m-2 s-1)':Gs_shade,
            'Pn (micromol plant-1 s-1)':Pn,  
            'Lp mean (kg m-2 s-1 kPa-1)':LpMean,
            'totalFlux (kg m-2 s-1)': flux_tot
            }
    
    T=np.zeros(ndays)
    LAI_out_day=np.zeros(ndays)
    Time_day=[0]*ndays
    i=0
    for d in drange:
        dend=int(d+24/timeRecord)
        T[i]= np.sum(Ep[d:dend])*timeRecord*3600
        LAI_out_day[i]=LAI[d]
        Time_day[i]=meteo.loc[d,'Time']
        i=i+1
        
    output_day={'Time':Time_day,
                'Es_wet (mm day-1)':Eswet,
                'Es_dry (mm day-1)':Esdry,
                'ETo (mm day-1)':ETo,
                'T (mm day-1)':T,
                'Irrigation (mm day-1)':IrriVol,
                'Rain (mm day-1)':RainDay,
                'Qd':Qd,
                'LAI':LAI_out_day
                }
    
    
    output=pd.DataFrame(output)
    
    output_day=pd.DataFrame(output_day)
    
    SWC_wet=pd.DataFrame(SWCwet)
    SWC_dry=pd.DataFrame(SWCdry)
    SWC=pd.concat([SWC_wet,SWC_dry],axis=1)
    column=meteo['Time']
    SWC.insert(loc=0,column='Time',value=column)
    
    dir_simulation_out=dir_simulation+'/output/'
    if not os.path.exists(dir_simulation_out):
        os.makedirs(dir_simulation_out)
        
    output.to_csv(dir_simulation_out+'out_instant.csv',index=False)
    output_day.to_csv(dir_simulation_out+'out_day.csv',index=False)
    SWC.to_csv(dir_simulation_out+'SWC.csv',index=False)
    
    final_time=datetime.datetime.now()
    runTime=final_time-ini_time   
    print('Elapsed  time: '+str(runTime))
    print('Program finished')
    
if __name__ == '__main__':  main()
