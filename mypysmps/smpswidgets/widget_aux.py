#!/usr/bin/env python
# -*- coding: utf-8 -*-
#################
import warnings
import sys
import os
from mypysmps.io.read import read

from IPython.display import Image, display, Math, Latex, HTML, Javascript
from ipywidgets import widgets, interact, Layout, interactive
#################
"""
mypysmps.smpswidgets.controlPanel
================

Auxiliary functions for interactive plots:
    checkfile

Created on Tue Jul 29 08:50 2020

@author: flovan / fvanden

Revision history:   29.07.2020 - Created


"""
## -------------------------------------------------------------------------- ##
##                                                                            ##
##                            Widget States                                   ##
##                                                                            ##
## -------------------------------------------------------------------------- ##

# mean free path
meanFreePath = widgets.BoundedFloatText( value=67.3,
        min=65.0,
        max=69.2,
        step=0.01,
        description='mean free path $nm$:',
        style = {'description_width': 'initial'},
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.2f',
        tooltip = ('Mean free path'))

# diameter
Dp = widgets.BoundedFloatText( value=14.3,
        min=4.37,
        max=750,
        step=0.1,
        description='Particle diameter $nm$:',
        style = {'description_width': 'initial'},
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',
        tooltip = ('Particle diameter'))

# Stokes number
Stk50 = widgets.BoundedFloatText( value=0.23,
    min=0,
    max=1.0,
    step=0.1,
    description='Stokes number:',
    style = {'description_width': 'initial'},
    disabled=False,
    continuous_update=True,
    orientation='horizontal',
    readout=True,
    readout_format='.2f',
    tooltip = ('Stokes number'))

# Particle Density
RhoP = widgets.BoundedFloatText( value=1.0000,
    min=0,
    max=1.5,
    step=0.1,
    description='Particle Density $g/cm3$:',
    style = {'description_width': 'initial'},
    disabled=False,
    continuous_update=True,
    orientation='horizontal',
    readout=True,
    readout_format='.1f',
    tooltip = ('Particle density'))

# Aerosol volumetric flow rate
Q = widgets.BoundedFloatText( value=0.373,
    min=0,
    max=1.,
    step=.001,
    description='Aerosol Flow Rate $L/min$:',
    style = {'description_width': 'initial'},
    disabled=False,
    continuous_update=True,
    orientation='horizontal',
    readout=True,
    readout_format='.3f',
    tooltip = ('Aerosol volumetric flow rate'))

# Sheath volumetric flow rate
Qc = widgets.BoundedFloatText( value=3.,
    min=0,
    max=10.,
    step=0.1,
    description='Sheath Flow Rate $L/min$:',
    style = {'description_width': 'initial'},
    disabled=False,
    continuous_update=True,
    orientation='horizontal',
    readout=True,
    readout_format='.3f',
    tooltip = ('Sheath volumetric flow rate'))

# gas viscosity
Eta = widgets.BoundedFloatText( value=1.83245*10**-4,
    min=1.0*10**-4,
    max=2.0*10**-4,
    step=1.*10**-6,
    description='Gas viscosity $g/cm*s$:',
    style = {'description_width': 'initial'},
    disabled=False,
    continuous_update=True,
    orientation='horizontal',
    readout=True,
    readout_format='.8f',
    tooltip = ('Reference gas viscosity'))

# diameter aerosol inlet
W = widgets.BoundedFloatText( value=0.0508,
    min=0.05,
    max=0.06,
    step=0.1,
    description='Nozzle diameter $cm$:',
    style = {'description_width': 'initial'},
    disabled=False,
    continuous_update=True,
    orientation='horizontal',
    readout=True,
    readout_format='.4f',
    tooltip = ('Diameter of aerosol inlet'))

# number of elementary charges particle
n = widgets.IntSlider( value=1,
    min=-6,
    max=6,
    step=1,
    description='E charges in prtcle $#$:',
    style = {'description_width': 'initial'},
    disabled=False,
    continuous_update=True,
    orientation='horizontal',
    readout=True,
    readout_format='.4f',
    tooltip = ('Number of elementary charges on particle'))

# elementary charge value
e = widgets.BoundedFloatText( value=1.6022*10**-19,
    min=1.*10**-19,
    max=2.*10**-19,
    step=.1*10**-19,
    description='elementary charge $A.s$:',
    style = {'description_width': 'initial'},
    disabled=False,
    continuous_update=True,
    orientation='horizontal',
    readout=True,
    readout_format='.4f',
    tooltip = ('Value of elementary charge'))

# viscosity select
etaselect = widgets.ToggleButtons(
    options=['Calculated', 'Reference'],
    description='Gas viscosity:',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    style =  {'description_width': 'initial'},
    tooltips=['Dynamic gas viscosity','Reference gas viscosity']
)

# free path select
freepathselect = widgets.ToggleButtons(
    options=['Calculated', 'Reference'],
    description='Mean free path:',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    style =  {'description_width': 'initial'},
    tooltips=['Dynamic mean free path','Reference mean free path']
)

# DMA select
DMAselect = widgets.ToggleButtons(
    options=['long', 'nano'],
    description='DMA:',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    style =  {'description_width': 'initial'},
    tooltips=['Long DMA','Nano DMA']
)

# pressure
Pressure = widgets.BoundedFloatText(
    value=101.3,
    min=100,
    max=110.0,
    step=0.1,
    description='Sample Pressure $kPa$:',
    style =  {'description_width': 'initial'},
    disabled=False
)

# temperature
Temperature = widgets.BoundedFloatText(
    value=296.15,
    min=200.0,
    max=300.0,
    step=1,
    description='Sample Temperature $K$:',
    style =  {'description_width': 'initial'},
    disabled=False
)


## -------------------------------------------------------------------------- ##
##                                                                            ##
##                            Widget Functions                                ##
##                                                                            ##
## -------------------------------------------------------------------------- ##



def wchooseFile(directory = "data/"):
    """
    """
    options = [os.path.join(root, name)
             for root, dirs, files in os.walk(directory)
             for name in files
             if name.endswith((".txt", ".csv"))]
    return sorted(options)

def wread(filename):
    """
    """
    fname, fextension = os.path.splitext(filename)
    
    if fextension == '.txt':
        fileorg = 'AIM_text'
    else:
        fileorg = 'AIM'
        
    SMPS = read(filename, fileorg)
    
    return SMPS
    
    
    


    