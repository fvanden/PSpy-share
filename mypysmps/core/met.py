# -*- coding: utf-8 -*-
#################
import datetime as dt

from ..config import get_metadata
from ..util.timetransform import TimeTransform
tt = TimeTransform()
#################

"""
mysmps.core.met
================

A general meteorological data class
    MET

Created on Wed Sep 30 12:17 2020

@author: flovan / fvanden

Revision history:   30.09.2020 - Created


"""
## -------------------------------------------------------------------------- ##


class MET(object):
    """
    A class storing meteorological data
    
    Parameters - required
    ---------------------
    time : dict
        Time of each sample
    date: dict
        Date of each sample
    data : dict
        Data fields of measurements
    
    Parameters - optional
    ---------------------
    
    
    """
    ## ------------------------------------------------------------------ ##
    ## Constructors/Destructors                                           ##
    ## ------------------------------------------------------------------ ##
    def __init__(self, **kwargs):

        self.__dict__.update(kwargs)
        self.instrument = 'MET'
    
    def __del__(self):
        pass


    ## ------------------------------------------------------------------ ##
    ## Methods                                                            ##
    ## ------------------------------------------------------------------ ##

    # private:
    
    # sample
    # header
    # metadata