# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 12:45:20 2019

@author: ogarcia
"""
import numpy as np
def sun_riseSet(doy,lat,long):
    """
    Fucntion to obtain sun rise and sunset hour at each time of the day and at
    each latitude. From Villalobos in ModelSPA V2.
    
    Inputs
    ------
    doy: day of year
    lat: latitude
    long: deviation in degrees from the standard meridian
    """
    latg=lat
    longg=long
    
    d=doy;
    latr = latg * np.pi / 180 #latitude in radians
    d1 = 2 * np.pi * (d - 1) / 365
    d2 = 2 * d1
    d3 = 3 * d1
    decr = 0.006918 - 0.399912 * np.cos(d1) + 0.070257 *np.sin(d1) - 0.006758 * np.cos(d2)+ 0.000907 * np.sin(d2) - 0.002697 * np.cos(d3) + 0.00148 * np.sin(d3) # sun declination (radian)
    e = 9.87 * np.sin(4 * (d - 81) / 364 * np.pi) - 7.53 * np.cos(2 * np.pi / 364 * (d - 81)) - 1.5 * np.sin(2 * np.pi / 364 * (d - 81))
    K1 = -np.tan(latr) * np.tan(decr);
    if K1 > 1:
        K1 = 1
    elif K1 < -1:
        K1 = -1
    elif K1 != 0:
         HSM1R = np.arctan(np.sqrt((1 - K1** 2) / K1** 2)) 
    else:
        HSM1R = 3.1416 / 2;
    
    if K1 < 0: 
        HSM1R = 3.1416 - HSM1R
    
    HSM1G = HSM1R * 180 / np.pi
    corre = 1.5 * 0.5 / np.cos(latr) * 24 / 360
    hsm = HSM1G / 7.5 + 2 * corre
    
    
    tsunr = 12 - hsm / 2 + longg / 15 - e / 60;
    tsuns = tsunr + hsm;
    
    return tsunr,tsuns 
#if __name__ == "__main__":
#    tsunr,tsuns=sun_riseSet(350,-34.95,4)