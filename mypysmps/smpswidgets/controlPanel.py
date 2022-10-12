#!/usr/bin/env python
# -*- coding: utf-8 -*-
#################
from IPython.display import Image, display, Math, Latex, HTML, Javascript
from ipywidgets import widgets, interact, Layout, interactive
import os
import numpy as np
import math
import matplotlib.pyplot as plt

from mypysmps.io.read import read
from mypysmps.graph.psdisplay import PSDisplay
from mypysmps.util.basic import D50, calcC, calcKn, conversion, particleFractions

#################
"""
mypysmps.smpswidgets.controlPanel
================

Control panel widget and functions for interactive plots:
    cpwidget

Created on Tue Jul 28 17:22 2020

@author: flovan / fvanden

Revision history:   28.07.2020 - Created


"""
## -------------------------------------------------------------------------- ##

def quicklook(SMPS):
    """
    Creates an interactive quicklook graph
    for SMPS data
    
    Parameters
    ----------
    SMPS : PSpy.core.smps object
        object containing SMPS data
        
    Returns
    -------
    ui3 : widget.VBox display
        widget that can be opened using IPython.display.display
        
    """
    # create particle sizer display class
    psd = PSDisplay(SMPS)
    
    # create sample slider
    sample = widgets.IntSlider(min=0, max=len(SMPS.sample['data'])-1, continuous_update = False, description = 'sample #',layout = Layout(width = '100%'))
    
    # create fields slider containing SMPS data variables
    sfields = widgets.ToggleButtons(options=SMPS.data['variables'], value = SMPS.data['variables'][0], description='Variable:', disabled=False, button_style='')
    
    # create D50 checkbox
    D50line = widgets.Checkbox(value=False,description='D50',disabled=False,indent=False, layout = Layout(width='20%'))
    
    # create boundary line checkbox
    boundaryline = widgets.Checkbox(value=False,description='Boundaries',disabled=False, indent=False, layout = Layout(width='20%'))
    
    # create periodical indicator checkbox
    indperiods = widgets.Checkbox(value=False,description='24 hour periods',disabled=False,indent=False, layout = Layout(width='20%'))
    
    # create starttime for periodical indicator text box
    starttime = widgets.Text(value='00:00:00', placeholder='HH:MM:SS', description='Start time:', disabled=False)
    
    # create count limit box
    clim = widgets.BoundedFloatText( value=2*10**4, min=1*10**1, max=1*10**8,step=1*10**1, description='clim:',
        style = {'description_width': 'initial'}, disabled=False, continuous_update=True, orientation='horizontal',
        readout=True, readout_format='.2f', layout = Layout(width='25%'))
    
    # create diameter limit box
    ylim = widgets.BoundedFloatText( value=1*10**3, min=1*10**1, max=1*10**8,step=1*10**1, description='dlim:', disabled=False, continuous_update=True, orientation='horizontal', readout=True, readout_format='.2f', layout = Layout(width='25%'))
    
    # create histogram outplot
    outplot =  widgets.interactive_output(psd.histogram, {'field':sfields,'sample':sample,'add_D50':D50line,'add_boundaries':boundaryline, 'xlim': ylim,'ylim':clim}) 
    
    # create heatplot outplot
    outplot2 = widgets.interactive_output(psd.plot, {'field':sfields,'clim':clim, 'ylim':ylim, 'indicator':sample,'periods':indperiods, 'starttime':starttime})
    
    # create printed date output
    date = widgets.interactive_output(psd.dateInfo,{'sample':sample})
    
    # create organised display widgets
    ui1 = widgets.HBox([sfields])
    ui2 = widgets.HBox([D50line, boundaryline, clim, ylim])
    uidate = widgets.HBox([date])
    uiperiods = widgets.HBox([indperiods, starttime])
    ui3 = widgets.VBox([ui1, ui2, outplot, sample, uidate, outplot2, uiperiods], layout = Layout(height='100%'))
    
    return ui3

def D50widget():
    """
    """
    

    

def show(lm, Dp, tprint = True):
    Kncalc = calcKn(lm, Dp)
    Ccalc = calcC(Kncalc)
    if tprint:
        print( ('Knudsen number: ' + '%.9f')%(Kncalc) )
        print( ('Cunningham Slip Correction: ' + '%.9f')%(Ccalc) )
    return Kncalc, Ccalc
    

def showFlow(Q):
    conv = conversion(Q,'L/min','cm3/s')
    print( ('Flow Rate (cm\N{SUPERSCRIPT THREE}/s): ' + '%.9f')%(conv) )

    
def showD50(lm, Dp, Stk50, RhoP, Q, Eta, W):
    Kncalc = calcKn(lm, Dp)
    C = calcC(Kncalc)
    Q = conversion(Q,'L/min','cm3/s')
    calcD50 = D50(Stk50, RhoP, Q, C, Eta, W)
    calcD50 = calcD50 * 10**7
    print( ('D50 (nm): ' + '%.9f')%(calcD50) )


def showZp(P,T,lm,Dp,n,Eta,fps,etas,e):
    if fps == 'Calculated':
        lm = freePathCalc(P, T, lambdar = lm, Pr = 101.3 , Tr = 296.15, S = 110.4)
    else:
        lm = lm
    if etas == 'Calculated':
        Eta = etaCalc(T, Tr = 296.15, S = 110.4, nr = conversion(Eta, 'g/cm-s','kg/m-s'))
    else:
        eta = Eta
        
    Kncalc = calcKn(lm, Dp)
    C = calcC(Kncalc)
    
    Zp = electricalMobility(n, C, Eta, Dp, e = 1.6022*10**-19)
    
    print( ('Zp E15 :'  + '%.9f')%(Zp*10**15) )
    
def showDMA(DMAselect):
    dmatype = _INSTRUMENT_SETTINGS['SMPS'][DMAselect]
    for key, value in dmatype.items():
        print(key,':', value)
    return None
     
def showVoltage(DMAselect, Eta, Qc,n,e, Dp, lm):
    _,C = show(lm, Dp, tprint=False)
    dmatype = _INSTRUMENT_SETTINGS['SMPS'][DMAselect]
    #L = conversion(dmatype['length'],'cm','nm')
    #r1 = conversion(dmatype['r1'], 'cm','nm')
    #r2 = conversion(dmatype['r2'],'cm','nm')
    L = dmatype['length']
    r1 = dmatype['r1']
    r2 = dmatype['r2']
    Qc = conversion(Qc, 'L/min','cm3/s')
    
    V = voltageCalc(L, r1, r2, Eta, Qc, n, e, Dp, C)
    print( ('Voltage:'  + '%.9f')%(V) )
    return V


    button = widgets.Button(description="Reset",tooltip = 'Reset to default values')
    button.on_click(widgets.interactive_output(run_all,{'lm': lm, 'Dp': Dp} ))

    outsw =  widgets.interactive_output(show,{'lm': lm, 'Dp': Dp} )
    outsx =  widgets.interactive_output(showD50,{'lm':lm,'Dp':Dp,'Stk50':Stk50,'RhoP':RhoP,'Q':Q,'Eta':Eta,'W':W} )
    conva = widgets.interactive_output(showFlow, {'Q':Q})
    convs = widgets.interactive_output(showFlow, {'Q':Qc})
    outsZp = widgets.interactive_output(showZp, {'P': Pressure, 'T':Temperature, 'lm':lm, 'Dp':Dp, 'n':n, 'Eta':Eta, 'fps':freepathselect, 'etas':etaselect, 'e':e})
    outsDMA = widgets.interactive_output(showDMA, {'DMAselect':DMAselect}) 
    outsV = widgets.interactive_output(showVoltage, {'DMAselect':DMAselect,'Eta':Eta,'Qc':Qc,'n':n,'e':e,'Dp':Dp,'lm':lm})

    ui = widgets.VBox([button, lm, Dp, Stk50, RhoP, Q, conva, Qc, convs, Eta, W ])
    ui2 = widgets.VBox([])
    ui3 = widgets.VBox([outsw, outsx, n, e, etaselect,freepathselect, Pressure, outsZp, DMAselect, outsDMA,outsV])
    uif = widgets.HBox([ui,ui2,ui3])

    display(uif)
    
