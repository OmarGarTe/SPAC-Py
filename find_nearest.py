# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 13:26:22 2019

@author: ogarcia
"""

import numpy as np
def find_nearest(array, value):
    """
    Find nearest value in a numpy array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]