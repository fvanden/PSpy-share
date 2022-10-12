# -*- coding: utf-8 -*-
#################
import numpy as np
#################

"""
mysmps.util.basic
================

Basic equations for SMPS particle sizer
    D50
    calcC
    calcKn
    particleFractions
    conversion
    etaCalc
    electricalMobility
    freePathCalc

Created on Tue Jul 22 10:01 2020

@author: flovan / fvanden

Revision history:   22.07.2020 - Created


"""
## -------------------------------------------------------------------------- ##

def D50(Stk50, RhoP, Q, C, Eta, W):
    """
    Returns the cut point diameter D50
    
    Parameters
    ----------
    Stk50 : float
        Stokes number for 50% collection efficiency = 0.23
        
    RhoP : float
        particle density (g/cm3)
        
    Q : float
        volumetric flow rate (cm3/s)
    
    C : float
        Cunningham Slip Correction
    
    Eta : float
        gas viscosity (g/(cm·s))
        
    W : float
        nozzle diameter (cm)
        
    Returns
    -------
    D50 : float
        cut point diameter
    """
    D50 = np.sqrt( (9*np.pi*Stk50*Eta*W**3)/(4*RhoP*C*Q))
    
    return D50

def calcC(Kn, a = 1.165, b = 0.483, c = 0.997):
    """
    Returns the Cunningham Slip Correction C
    
    Parameters
    ----------
    Kn : float
        Knudsen Number (Kim et al., 2005)
        
    a : float
        fixed value
        
    b : float
        fixed value
    
    c : float
        fixed value
        
    Returns
    -------
    C : float
        Cunningham Slip Correction
    
    """
    C = 1 + Kn*( a + b * np.exp( - (c/Kn) )  )
    
    return C

def calcKn(lm, Dp):
    """
    Returns the Knudsen Number (Kim et al., 2005)
    
    Parameters
    ----------
    lm : float
        gas mean free path (nm)
        
    Dp : paticle diameter (nm)
        
    Returns
    -------
    Kn : float
        Knudsen Number
    """
    #lm = lm *10**-8
    #Dp = Dp *10**-8
    
    Kn = (2*lm) / Dp
    
    return Kn

def particleFractions(Dp, N):
    """
    """
    coeffDict = {-2: {0:-26.3328,1:35.9044,2:-21.4608,3:7.0867,4:-1.3088,5:0.1051},
                 -1: {0:-2.3197,1:0.6175,2:0.6201,3:-0.1105,4:-0.1260,5:0.0297},
                  0: {0:-0.0003,1:-0.1014,2:0.3073,3:-0.3372,4:0.1023,5:-0.0105},
                  1: {0:-2.3484,1:0.6044,2:0.4800,3:0.0013,4:-0.1553,5:0.0320},
                  2: {0:-44.4756,1:79.3772,2:-62.8900,3:26.4492,4:-5.7480,5:0.5049}
    }
    
    outs = np.empty([len(Dp), 5])
    outs[:,:] = np.nan
    
    for i in range(0, 5):
        coeff = coeffDict[N][i]
        outs[:,i] = (coeff * np.log(Dp)**i)
        
    fout = 10*np.sum(outs, axis = 1)    
    
    return fout

def etaCalc(T, Tr = 296.15, S = 110.4, nr = 1.83245*10**-5):
    """
    Calculates dynamic gas viscosity in kg*m-1*s-1
    
    Parameters
    ----------
    T : float
        Temperature (K)
        
    Tr : float
        Reference Temperature (K)
    
    S : float
        Sutherland constant (K)
        
    nr : float
        Reference dynamic viscosity
        
    Returns
    -------
    eta : float
        Dynamic gas viscosity in kg*m-1*s-1
        
    """
    eta = nr * ( (Tr + S) / (T+S) )*(T/Tr)**(3/2)
    return eta

def freePathCalc(P, T, lambdar = 6.730*10**-8, Pr = 101.3 , Tr = 296.15, S = 110.4):
    """
    Calculates gas mean free path in m
    
    Parameters
    ----------
    P : float
        pressure
        
    T : float
        Temperature (K)
        
    lambdar : float
        reference mean free path (m)
        
    Pr : float
        Reference pressure (kPa)
    
    Tr : float
        Reference Temperature (K)
    
    S : float
        Sutherland constant (K)
        
    Returns
    -------
    lambdao : float
        Gas mean free path in m
    """
    
    lambdao = lambdar * (Pr/P)*(T/Tr)*((1+S/Tr)/(1+S/T))
    
    return lambdao

def electricalMobility(n, C, eta, Dp, e = 1.6022*10**-19):
    """
    Calculates dynamic gas viscosity in kg*m-1*s-1
    
    Parameters
    ----------
    n : float
        number of elementary charges on a particle
        
    e : float
        elementary charge (1.6022 x 10-19 As)
    
    C : float
        Cunningham slip correction 
        
    eta : float
        dynamic gas viscosity
    
    Dp : float
        particle diameter
        
    Returns
    -------
    Zp : float
        electrical mobility
    """
    Zp = (n*e*C)/(3*np.pi*eta*Dp)
    
    return Zp

def voltageCalc(L, r1, r2, Eta, qc, n, e, Dp, C):
    """
    """
    IC = L / np.log(r1/r2)
    V = ((3 * Eta * qc)/(2*n*e*IC)) * (Dp/C)
    return V

def conversion(measure, fromm, tom):
    """
    Converts values from one units to another
    
    Parameters
    ----------
    measure : float
        value to be converted
        
    fromm : str
        units of value to be converted
    
    tom : str
        units to convert to 
        
    Returns
    -------
    outm : float
        converted value
    """
    # 1 L/min == 16.67 cm3/s, 1000/60
    if fromm == 'L/min' and tom == 'cm3/s':
        outm = measure *1000/60
    if fromm == 'L/min' and tom == 'm3/s':
        outm = measure/1000/60
       
    # gas viscosity
    if fromm == 'Pa*s' and tom == 'g/cm-s':
        outm = measure * 10
    if fromm == 'g/cm-s' and tom == 'Pa*s':
        outm = measure / 10
    if fromm == 'kg/m-s' and tom == 'g/cm-s':
        outm = measure * 10
    if fromm == 'g/cm-s' and tom == 'kg/m-s':
        outm = measure / 10
    if fromm == 'kg/m-s' and tom == 'Pa*s':
        outm = measure
    if fromm == 'Pa*s' and tom == 'kg/m-s':
        outm = measure
        
    # temperature
    if fromm == 'celcius' and tom == 'kelvin':
        outm = measure + 273.15
    if fromm == 'kelvin' and tom == 'celcius':
        outm = measure - 273.15
        
    # distance
    if fromm == 'nm' and tom == 'm':
        outm = measure * (10**-9) 
    if fromm == 'm' and tom == 'nm':
        outm = measure * (10**9) 
    if fromm == 'nm' and tom == 'cm':
        outm = measure / (10**-7)   
    if fromm == 'cm' and tom == 'nm':
        outm = measure * (10**7)     
    if fromm == 'nm' and tom == 'microm':
        outm = measure / 1000    
    if fromm == 'microm' and tom == 'nm':
        outm = measure * 1000 
    if fromm == 'cm' and tom == 'm':
        outm = measure / 100
    if fromm == 'm' and tom == 'cm':
        outm = measure * 100
        
    return outm
    

    