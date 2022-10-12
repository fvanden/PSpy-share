#!/usr/bin/env python
# -*- coding: utf-8 -*-
#################
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from ipywidgets import widgets, GridspecLayout, Button, Layout, interact, interactive
from IPython.display import Image, display, Math, Latex, HTML, Javascript

#from mypysmps.io.read import read
#################
"""
mypysmps.smpswidgets.app_aux
================

Auxiliary functions for quicklook app:


Created on Thu Feb 4 14:50 2021

@author: flovan / fvanden

Revision history:   04.02.2021 - Created


"""
## -------------------------------------------------------------------------- ##
##                                                                            ##
##                                  Functions                                 ##
##                                                                            ##
## -------------------------------------------------------------------------- ##


def getDefaultValues(IDict,instr_select):
    """
    gets Default values
    """
    instrument = IDict[instr_select]
    if instrument is None:
        variable_list =['var1','var2','var3']
        default_value = 'var1'
        diameter_list = [0,1,2,3]
        diamdefval = 0
        maxSample = 100
        fileorg = 'none'
    else:

        variable_list = np.append(instrument.data['variables'],list(instrument.__dict__.keys()))
        maxSample = len(instrument.sample['data'])-1

        diameter_list = instrument.diameter['data']
        diamdefval = diameter_list[0]

        if instrument.instrument_type == 'OPC':
            default_value = 'raw_counts'
            fileorg = 'OPC'
        elif instrument.instrument_type == 'SMPS':
            default_value = 'total_concentration'
            fileorg = 'AIM'
        else:
            default_value = 'total_concentration'
            fileorg = 'AIM'

    return fileorg, default_value, variable_list, diameter_list, diamdefval, maxSample


def update_var_widget(IDict,instr_select):
    """
    """

    fileorg, default_value, variable_list,maxSample = getDefaultValues(IDict, instr_select)

    variable.set_trait('options', variable_list)
    variable.set_trait('value',default_value)

    sample.set_trait('max',maxSample)

    return variable
