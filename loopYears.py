# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 11:16:09 2021

@author: omarg

Loop through different years
This is a script dveloped to compare irrigation strategies using past or forecast
weather data. 

The variable weather calls the weather file that contains all the years in a row
The variable irri call the irrigation file. If the irrigation file is based
on past data it is named: "Trad_IrriFile_model.csv", on the contrary, if 
the irri file is based on forecast data it is named: "WRF_IrriFile.csv"

The script save the results for each year in the corresponding output file, and
the results for all the years toughether. 
"""

import pandas as pd
import numpy as np
import main_energy


years=np.arange(2013,2020)

RDI_opt=False#reduce irrigation volume during a desired period
iniRDI=196#initial doy for RDI
endRDI=258#end doy for RDI
K_RDI=0.3#reduction coeficient

file=open('experiment_folder_name.txt','r')
exp_folder=str(file.read())
file.close()



#Files with various years

weather=pd.read_csv(exp_folder+'weatherAlmond_model_Rad_modify.csv')
irri=pd.read_csv(exp_folder+'WRF_corr_IrriFile_model.csv')

weather['Time']=pd.to_datetime(weather['Time'])
weather['Year']=weather['Time'].dt.year


irri['Time']=pd.to_datetime(irri['Time'])
irri['Year']=irri['Time'].dt.year





#Read files to be updated
bb_date=pd.read_csv(exp_folder+'/input/tree/budbreakdates.csv')
lai=pd.read_csv(exp_folder+'/input/tree/LAI.csv')
rld=pd.read_csv(exp_folder+'/input/tree/rld.csv')

#Container for the data
out_inst_all=pd.DataFrame([])
out_day_all=pd.DataFrame([])

for y in years:
    
    #Slice over files
    
    weatherY=weather[weather['Year']==y]
    weatherY.reset_index(inplace=True,drop=True)

    
    #Slice irri file
    irriY=irri[irri['Year']==y]
    irriY.reset_index(inplace=True,drop=True)
    
    if RDI_opt:
        irriY['DOY']=irriY.loc[:,'Time'].dt.dayofyear
        irriY.set_index('DOY',inplace=True)
        irriY.loc[iniRDI:endRDI,'Irrivol']=irriY.loc[iniRDI:endRDI,'Irrivol']*K_RDI
        irriY.reset_index(inplace=True,drop=True)
        
    
    #Update files dates
    
    weatherY.to_csv(exp_folder+'/input/meteo/weather.csv',index=False)
    irriY.to_csv(exp_folder+'/input/irrigation/IrriFile.csv',index=False)
    
    bb_date['Date']=weatherY.loc[0,'Time']
    bb_date.to_csv(exp_folder+'/input/tree/budbreakdates.csv',index=False)
    
    lai.loc[0,'Date']=weatherY.loc[0,'Time']
    lastDay=np.max(weatherY['DOY'])
    rowLastDay=np.where(weatherY['DOY']==lastDay)
    lai.loc[1,'Date']=weatherY.loc[rowLastDay[0][0],'Time']
    lai.to_csv(exp_folder+'/input/tree/LAI.csv',index=False)
    
    rld.loc[:,'Time']=weatherY.loc[0,'Time']
    rld.to_csv(exp_folder+'/input/tree/rld.csv',index=False)
    
    #Run the model
    main_energy.main()
    print(y)
    #Create outputs per year
    out_inst=pd.read_csv(exp_folder+'\output\out_instant.csv')
    out_inst.to_csv(exp_folder+'\output\out_instant_'+str(y)+'.csv',index=False)
    out_inst_all=out_inst_all.append(out_inst)
    out_day=pd.read_csv(exp_folder+'\output\out_day.csv')
    out_day.to_csv(exp_folder+'\output\out_day_'+str(y)+'.csv',index=False)
    out_day_all=out_day_all.append(out_day)

out_inst_all.to_csv(exp_folder+'\output\out_instant_all.csv',index=False)
out_day_all.to_csv(exp_folder+'\output\out_day_all.csv',index=False)
