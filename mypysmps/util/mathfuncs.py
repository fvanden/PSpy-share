# -*- coding: utf-8 -*-
#################
import numpy as np
#################

"""
mysmps.util.mathfuncs
================

Basic mathematical equations
    roundup

Created on Thu Aug 6 18:00 2020

@author: flovan / fvanden

Revision history:   06.08.2020 - Created


"""
## -------------------------------------------------------------------------- ##

def roundup(x):
    """
    Rounds value up to nearest 10, 100, 1000, ...
    
    Parameters
    ----------
    x : float
        input number
        
    Returns
    -------
    y : float
        x rounded up to nearest 10, 100, 1000 etc
    """
    d = np.round(100*np.log(x)/np.log(10))/100
    d = np.floor(d)
    y = int(np.ceil(x / 10**d)) * 10**d
    return(y)