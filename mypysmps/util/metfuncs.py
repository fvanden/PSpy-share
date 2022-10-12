# -*- coding: utf-8 -*-
#################
import numpy as np
import math
#################

"""
mysmps.util.metfuncs
==================

Functions for working with MET data


Created on Mon Nov 29 14:04 2021

@author: flovan / fvanden

Revision history:   29.11.2021 - Created

"""


def averageWind(windspeed, winddirection, out = 'optimal'):
    """
    Average all windspeed and winddirection values in given
    lists

    Parameters
    ----------
    windspeed : list of floats
        windspeed values in m/s

    winddirection : list of floats
        wind direction in degrees

    optimal : bool
        if set to True, the scalar average is calculated
        for the wind speed and the vector average is used
        for the wind directions

    Returns
    -------
    aws : float
        average wind speed for given list of values

    awd : float
        average wind direction for given list of values
    """

    u,v = calcScalarWind(windspeed, winddirection)

    aws = np.average(windspeed)

    avws, avwd = averageWindSpeedDirection(u,v)

    if out == 'optimal':
        return aws, avwd
    else:
        return avws, avwd

def calcScalarWind(windspeed, winddirection):
    """
    Calculates vectorial wind from windspeed and
    directions

    Parameters
    ----------
    windspeed : list of floats
        windspeed in m/s

    winddirection : list of floats
        wind direction in degrees

    Returns
    -------
    u : list of floats
        u (E-W) component vectorial wind

    v : list of floats
        v (N-S) component vectorial wind

    """
    u = -np.asarray(windspeed) * np.sin( 2*np.pi * np.asarray(winddirection)/360 )
    v = -np.asarray(windspeed) * np.cos( 2*np.pi * np.asarray(winddirection)/360 )

    return u, v

def scalar2regular(u,v):
    """
    Calculates windspeed and direction in m/s?
    and degrees

    Parameters
    ----------
    u : list of floats
        u (E-W) component vectorial wind

    v : list of floats
        v (N-S) component vectorial wind

    Returns
    -------
    windspeed : list of floats
        windspeed in m/s

    winddirection : list of floats
        wind direction in degrees (North is always 0.0)
    """

    windspeed = np.sqrt(np.asarray(u)**2 + np.asarray(v)**2)
    winddirection = (270 - np.arctan2(np.asarray(v),np.asarray(u))*(180/np.pi))%360

    return windspeed, winddirection



def averageWindSpeedDirection(u,v):
    """
    Calculates average wind speed and direction from
    vectorial wind components

    Parameters
    ----------
    u : list of floats
        u (E-W) component vectorial wind

    v : list of floats
        v (N-S) component vectorial wind

    Returns
    -------
    avws : float
        average wind speed (m/s)

    avwd : float
        average wind direction (degrees)
    """
    avws = (np.average(u)**2 + np.average(v)**2)**0.5
    avwd = (math.atan2(np.average(u),np.average(v))*(360/2/np.pi)) + 180

    return avws, avwd

def calculateTrueWind(A,V,AWD,alpha):
    """
    Calculates the True Wind from apparent wind and ships speeds and directions

    Can be validated by entering values from for example this website:
    https://www.nauticed.org/freesailingcourse-m1-2

    Parameters
    ----------
    A : float
        Apparent wind speed (consistent units for all speeds)

    V : float
        Ship speed on ground (SOG) - (consistent units for all speeds)

    AWD : float
        Apparent wind direction (degrees):
        The direction from which the wind is blowing, as measured/observed on
        the ship (in degrees from True North)

    alpha : float
        Ship heading (degrees from True North)

    Returns
    -------
    T : float
        True wind speed (in given units)

    TWD : float
        True wind direction (degrees):
        The direction from which the wind is blowing, corrected for the ships
        heading and speed (in degrees from True North)

    """

    tango = (AWD - alpha)%360

    T = np.sqrt(V**2 + A**2 - 2*V*A*np.cos(np.deg2rad(tango))  )

    if 0 < tango < 180: # starboard
        TWD = np.rad2deg(np.arccos((A * np.cos(np.deg2rad(tango)) - V) / T))%360
    else: # port
        TWD = -np.rad2deg(np.arccos((A * np.cos(np.deg2rad(tango)) - V) / T))%360

    TWD = (TWD + alpha)%360

    return T, TWD
