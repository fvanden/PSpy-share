# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 10:49:59 2018

@author: fvanden


"""
from scipy import signal
from scipy import ndimage
from scipy import fftpack
import numpy as np

import cmath

from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

#from Fourier.misc import removeNans


################################

# FRACVAR SCRIPTS

################################


def createFreqPlot(ax, s1, s3, Fs1 = 25., Fs2 = 25., config = None, markersize = 20, fracCorr2 = 1):
    """
    Creates a stem plot with Fractions of variance 
    explained based on Fourier spectral analysis
    for two different series of data
    
    Parameters
    ----------
    ax : matplotlib figure axis
        axis of the figure on which to plot
        
    s1 : list of floats or ints
        data series (Fourier transform is performed on this)
        
    s3 : list of floats
        data series (Fourier transform is performed on this)
        
    Fs1 : float
        frequency of s1 data series
    
    Fs2 : float
        frequency of s3 data series
        
    config: dict
        dictionary with plot info
        
    markersize : int 
        size of the markers for the plot
    
    fracCorr2 : int 
        a correction factor ?
        
    """
    plt.axes(ax)
    
    if config is None:
        plot1lab = 'real ML'
        plot2lab = 'artificial ML'
        xlab = 'Distance lag (m)'
        rescor1 = 1
        rescor2 = 1
    else:
        plot1lab = config['plot1lab']
        plot2lab = config['plot2lab']
        if config['rescor1'] and config['rescor2']:
            rescor1 = config['rescor1']
            rescor2 = config['rescor2']
            units = config['units1'] + '_' + config['units2']
        elif config['rescor1'] == 60:
            units = 'hours'
            rescor1 = config['rescor1']
            rescor2 = config['rescor1']
        else:
            units = config['units']
            rescor1 = config['rescor1']
            rescor2 = config['rescor1']
        xlab = config['xlab'] + '[' + units + ']'
    
    freqs1, spds1,fraqvars1 = computeFourierInfo(s1, Fs1)
    freqs3, spds3,fraqvars3 = computeFourierInfo(s3, Fs2)
    
    fraqvars3 = np.asarray(fraqvars3) * fracCorr2
    
    freqs = np.matrix([ ((1/np.asarray(freqs3))/rescor2)[1:11], ((1/np.asarray(freqs1))/rescor1)[1:11]])
    

    ax.set_title('Fractions of explained variance ')
    
    markerline, stemlines, baseline = ax.stem(np.arange(10), fraqvars1[1:11], 'k--', label = plot1lab, markerfmt="ko")
    color = 'xkcd:green'
    plt.setp(stemlines, color = color)
    plt.setp(markerline, color = color, markeredgecolor = color, markersize = markersize)
    plt.setp(baseline, color='k')
    
    markerline, stemlines, baseline = ax.stem(np.arange(0.2,10), fraqvars3[1:11], 'k--', label = plot2lab, markerfmt="ko")
    color = 'xkcd:sky blue'
    plt.setp(stemlines, color = color)
    plt.setp(markerline, color = color, markeredgecolor = color, markersize = markersize)
    plt.setp(baseline, color='k')
    
    for i in np.arange(.5,10):
        plt.plot([i,i],[0., 0.3], 'k--')
    
    lab = [('%1.1f - \n %1.1f')%(x,y) for x, y in zip(freqs.max(axis=0).tolist()[0], freqs.min(axis=0).tolist()[0])]
    
    plt.xticks(np.arange(0.1,10),lab, rotation = 'vertical') 
    
    ax.set_xlim([-0.5,10])
    ax.set_ylim([0.,0.3])
    plt.ylabel(r'$|P(f)|$' + r'$^2$' + r'$/\sigma_{f}^{2}$')
    plt.xlabel(xlab)
    ax.legend(loc='upper center')
    
    return ax
    
def createSpectraPlot(ax, s1, s3, Fs = 25, config = None):
    """
    """
    
    if config is None:
        plot1lab = 'real ML'
        plot2lab = 'artificial ML'
        xlab = 'Distance lag (m)'

    else:
        plot1lab = config['plot1lab']
        plot2lab = config['plot2lab']
        xlab = config['xlab']
        units = config['units']


        
    def func(x,a,b):
        return a+x*b
        
    st = 0
    
    
    freq1, psd1,fraqvars1 = computeFourierInfo(s1, Fs)
    freq1 = np.asarray(freq1)[~np.isnan(np.asarray(freq1))]
    psd1 = np.asarray(psd1)[~np.isnan(np.asarray(psd1))]
    
    ed = int(np.round((0.12*len(freq1))))
    
    # check if freq values are always decreasing
    
    idxs = np.where(np.diff(np.asarray(1/freq1))>0)
    if len(idxs[0]) == 0:
        idx = -1
    elif len(idxs[0]) == 1:
        idx = idxs[0][0]
    else:
        print("folding frequencies")
        return
    
    freq1 = freq1[:idx]
    psd1 = psd1[:idx]

    ax.loglog(freq1, psd1, color = 'xkcd:green' ,label = plot1lab)
    popt1, pcov1 = curve_fit(func, np.log10(freq1[st:ed]), np.log10(psd1[st:ed]))
    ax.loglog(freq1[st:ed], 10**(np.log10(freq1[st:ed])*popt1[1] + popt1[0]), 'k--')
    sigma1 = np.sqrt([pcov1[0,0], pcov1[1,1],]) 
   
    freq3, psd3,fraqvars3 = computeFourierInfo(s3, Fs)
    freq3 = np.asarray(freq3)[~np.isnan(np.asarray(freq3))]
    psd3 = np.asarray(psd3)[~np.isnan(np.asarray(psd3))]
    
    # check if freq values are always decreasing
    
    idxs = np.where(np.diff(np.asarray(1/freq3))>0)
    if len(idxs[0]) == 0:
        idx = -1
    elif len(idxs[0]) == 1:
        idx = idxs[0][0]
    else:
        print("folding frequencies")
        return
    
    freq3 = freq3[:idx]
    psd3 = psd3[:idx]
    
    ax.loglog(freq3, psd3, color = 'xkcd:sky blue' ,label = plot2lab)
    popt3, pcov3 = curve_fit(func, np.log10(freq3[st:ed]), np.log10(psd3[st:ed]))
    ax.loglog(freq3[st:ed], 10**(np.log10(freq3[st:ed])*popt3[1] + popt3[0]), 'k--')
    sigma3 = np.sqrt([pcov3[0,0], pcov3[1,1],]) 
    
    ax.set_xlabel((r'$\log_{10}$' + 'f(1/%s)')%(units))
    ax.set_ylabel(r'$\log_{10}$' + '(|P(f)|' + r'$^2$' + ')')
    ax.legend(loc='upper center')
    
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    
    xlab = xlims[0] + 0.1*xlims[0]
    ylab1 = ylims[0] + 5*ylims[0]
    ylab2 = ylims[0] + 3*ylims[0]
   
    ax.text(xlab,ylab1, ((r'$ \beta $ : %.2f $\pm$ %.2f, %s')%(np.abs(popt1[1]),sigma1[1], plot1lab)) )
    ax.text(xlab,ylab2, ((r'$ \beta $ : %.2f $\pm$ %.2f, %s')%(np.abs(popt3[1]),sigma3[1], plot2lab)) )
    
    ax.set_title('Spectra')
    
    return ax


################################

# PREPARATION SCRIPTS

################################

def belltaper(x):
    xlen = len(x)
    wind = np.ones(xlen)
    wind[:int(0.1*xlen)] = np.sin(5 * np.pi * np.arange(int(0.1*xlen))/xlen)**2
    wind[int(0.9*xlen)-1:] = np.sin(5 * np.pi * np.arange(int(0.9*xlen+1)-1, int(xlen+1))/xlen)**2
    
    return x * wind
    
def removeNansSimple(baselist, xlist):
    newlist = []
    newx = []
    for i in range(0, len(baselist)):
        if np.isnan(baselist[i]):
            pass
        else:
            newlist.append(baselist[i])
            newx.append(xlist[i])
    return newlist, newx

def my_detrend_simple(mldata,xdata):
    """
    """
    def func(x,a,b):
        return a+x*b

    popt, pcov = curve_fit(func,xdata, mldata)
    lintrend = np.asarray(xdata) * popt[1] + popt[0]
    data_out = mldata - lintrend
    
    return data_out, lintrend, popt[1] 
    
def my_detrend(mldata, xdata, idx = None):
    """
    """
    if idx is None:
        idx = np.argmin(np.abs(mldata))

    def func(x,a,b):
        return a+x*b

    popt1, pcov1 = curve_fit(func,xdata[:idx], mldata[:idx])
    lintrend1 = np.asarray(xdata[:idx]) * popt1[1] + popt1[0]
    data_out1 = mldata[:idx] - lintrend1

    popt2, pcov2 = curve_fit(func,xdata[idx:], mldata[idx:])
    lintrend2 = np.asarray(xdata[idx:]) * popt2[1] + popt2[0]
    data_out2 = mldata[idx:] - lintrend2

    data_out3 = np.append(data_out1, data_out2)
    lintrend_out = np.append(lintrend1, lintrend2)
    
    return data_out3, lintrend_out


def prepareData(mldata, xdata, lim = 80, detrending = None):
    """
    Prepares and conditions data series for Fourier transform
    
    Inputs
    ------
    mldata, list of floats : the series
    
    xdata, list of floats : the x data of the series
    
    Outputs
    -------
    x: x data with Nans removed for locations where ml data are Nans
    
    ml: ml data with Nans removed
    
    tdml: belltapered, detrended ml data, Nans removed
    
    medml: median filtered ml data, Nans removed
    
    tdmedml: belltapered, median filtered, detrended ml data, Nans removed
    
    lintrend: identified linear trend
    
    slope: slope of the linear trend
    
    """
    ml, x = removeNans([mldata], [xdata])
    ml = ml[0]
    x = x[0]
    medml = ndimage.median_filter(ml,25)
    if len(ml) > lim :
        if detrending is None:
            tdml,lintrend,slope  = my_detrend_simple(ml, x)
        elif detrending == 1:
            tdml,lintrend  = my_detrend(ml, x, idx = 576)
            slope = None
            
        tdml = belltaper(tdml)    
        tdmedml,_,_ = my_detrend_simple(medml, x)
        tdmedml = belltaper(tdmedml)
    
        return x, ml, tdml, medml, tdmedml, lintrend, slope
    else:
        print(len(ml))
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

################################

#      FOURIER SCRIPTS

################################

# the Fourier transform

def computeF(x, Fs):
    """
    Fast Fourier transform
    """
    L = len(x)
    qnp = np.asarray(x, dtype = float)
    

    fP1NoShift = np.fft.fft(qnp/L)
    
    ffreqNoShift = fftpack.fftfreq(L, d=float(Fs))
    freq = ffreqNoShift[:int(L/2)]
    
    return ffreqNoShift,fP1NoShift

def myFourierS(A):
    """
    Fourier transform based on R. Stull
    """
    x = np.asarray(A, dtype = float)    
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N,1))
    
    Fa = []
    da = []
    
    for j in k: 
        ip = 0   
        for i in range(0, N):
            ip += (1/N) * A[i] * np.cos(2 * np.pi * i * j / N) - (1j / N) * A[i] * np.sin(2 * np.pi * i * j /N) 
        Fa.append((ip[0].real, ip[0].imag))
        da.append(ip[0])
        
    return Fa, da

def myInvFourier(F):
    """
    Inverse Fourier transform based on R.Stull
    """
    N = len(F)
    n = np.arange(N)
    k = n.reshape((N,1))
    
    ip = 0
    Ak = np.zeros(N)
    
    for j in k:
        ip = 0
        for i in range(0, N):
            ip += (F[i].real * np.cos(2 * np.pi * i * j / N) - F[i].imag * np.sin( 2 * np.pi * i * j /N))
        Ak[j] = ip        
        
    return Ak

def sqrnrm(Fa):
    """
    Calculates the square of the norm
    """
    sqrn = []
    for i in range(0,len(Fa)):
        if i == 0:
            sqrn.append(np.nan)
        else:
            sqrn.append( Fa[i].real**2 + Fa[i].imag**2 )
    return sqrn

def discSpectDens(Fa, freqs):
    """
    Calculates the discrete spectral density
    """
    spd = []
    f = []
    if len(Fa) % 2 == 0:
        for i in range(0,len(Fa)):
            if i == 0:
                spd.append(np.nan)
                f.append(np.nan)
            elif i == int(len(Fa)/2):
                spd.append((Fa[i].real**2 + Fa[i].imag**2 ))
                f.append(freqs[i])
            elif i <= int(len(Fa)/2):
                spd.append( 2 *(Fa[i].real**2 + Fa[i].imag**2 ))
                f.append(freqs[i])
            else:
                spd.append(np.nan)
                f.append(np.nan)
    else:
        for i in range(0, int(len(Fa))):
            if i == 0:
                spd.append(np.nan)
                f.append(np.nan)
            elif i <= int(len(Fa)/2):
                spd.append( 2 *(Fa[i].real**2 + Fa[i].imag**2 ))
                f.append(freqs[i])
            else:
                spd.append(np.nan)
                f.append(np.nan)
                
    return spd, f

def specDens(spd, dn):
    """
    """
    Sa = []
    for i in range(0, len(spd)):
        Sa.append(spd[i]/dn)
    return Sa

def computeFourierInfo(tds, F = 25., message = False):
    """
    """
    # run transform
    freq, da = computeF(tds, F)
    
    # calculate the square of the norm
    Fasq = sqrnrm(da)
    
    # calculate the discrete spectral density
    spd, freq = discSpectDens(da, freq)
    
    # calculate the spectral density
    sa = specDens(spd, 1.)
    
    # check if results are consistent
    if message is True:
        print('First component equals mean:', np.mean(tds) ,  da[0] )
        print('sum spectral density equals variance:', np.nansum(sa) , np.var(tds) )
    
    return freq, spd, Fasq/np.var(tds)
    
def reduceDim(x, d = 4):
    """
    """
    Fa, da = myFourierS(x)
    adda = np.asarray(da)
    adda[d+1:] = 0
    adda[-d:] = da[-d:]
    adda[0] = da[0]
    mysignal = myInvFourier(adda)
    
    return mysignal
    
def reduceDimFast(x, d = 4, Fs = 25):
    """
    """
    Fa, da = computeF(x, Fs = Fs)
    adda = np.asarray(da)
    adda[d+1:] = 0
    adda[-d:] = da[-d:]
    adda[0] = da[0]
    mysignal = np.fft.ifft(adda)
    #mysignal = mysignal + da[0]
    mysignal = mysignal * len(x)
    
    return mysignal
    
def replaceCompFast(x1,x2, d = 4, Fs = 25, m= True, fracCorr2 = 1):
    """
    """
    Fa1, da1 = computeF(x1, Fs = Fs)
    Fa2, da2 = computeF(x2, Fs = Fs)
    adda1 = np.asarray(da1)
    adda2 = np.asarray(da2)
    
    da2 = np.asarray(da2) * fracCorr2
    
    #adda1[d+1:] = 0
    adda1[-d:] = da2[-d:]
    if m is True:
        adda1[0] = da2[0]
    else:
        adda1[0] = da1[0]
    mysignal = np.fft.ifft(adda1)
    
    mysignal = mysignal * len(x1)
    
    return mysignal
    
def replaceSingleComp(x1,x2, d = [4], Fs = 25, m= True, fracCorr2 = 1):
    """
    """
    Fa1, da1 = computeF(x1, Fs = Fs)
    Fa2, da2 = computeF(x2, Fs = Fs)
    adda1 = np.asarray(da1)
    adda2 = np.asarray(da2)
    
    da2 = np.asarray(da2) * fracCorr2
    
    #adda1[d+1:] = 0
    for ds in d:
        adda1[ds] = da2[ds]
        
    if m is True:
        adda1[0] = da2[0]
    else:
        adda1[0] = da1[0]
    mysignal = np.fft.ifft(adda1)
    
    mysignal = mysignal * len(x1)
    
    return mysignal
    
def replaceSingleAmplitude(x1, amps, d=[4], Fs = 25, m = True):
    """
    """
    Fa1, da1 = computeF(x1, Fs = Fs)
    adda1 = np.asarray(da1)
    
    for ds in d:
        amp, ph = cmath.polar(adda1[ds])
        #print(amp, ph)
        adda1[ds] = cmath.rect(amps[ds], ph)
        #print(cmath.rect(amps[ds], ph))
        
    mysignal = np.fft.ifft(adda1)
    
    mysignal = mysignal * len(x1)
    
    return mysignal
    
def replaceAmplitudes(x1, amps, d=4, Fs = 25, m = True):
    """
    """
    Fa1, da1 = computeF(x1, Fs = Fs)
    adda1 = np.asarray(da1)
    
    if m is True:
        s = 0
    else:
        s = 1
    
    for i in range(s, len(adda1[:d])):
        ds = i
        #print(i)
        amp, ph = cmath.polar(adda1[ds])
        #print(amp, ph)
        adda1[ds] = cmath.rect(amps[ds], ph)
        #print(cmath.rect(amps[ds], ph))
                
    mysignal = np.fft.ifft(adda1)
    
    mysignal = mysignal * len(x1)
    
    return mysignal
    
def replaceAmplitudeAtFrequency(x1, amps, frecbins, d = [4], Fs = 25):
    """
    """
    Fa1, da1 = computeF(x1, Fs = Fs)
    #print(Fa1)
    adda1 = np.asarray(da1)
    
    for ds in d:
        frec = Fa1[ds]
        amp, ph = cmath.polar(adda1[ds])
        for i in range(0, len(frecbins)-1):
            #print(frec, frecbins[i])
            if frec <= frecbins[i+1] and frec > frecbins[i]:
                idxbin = i
                break
            else:
                idxbin = None
                
        
        if idxbin is None:
            pass
        else:
            print(frec, idxbin, amps[idxbin],frecbins[idxbin])
            adda1[ds] = cmath.rect(amps[idxbin], ph)
        
    mysignal = np.fft.ifft(adda1)
    
    mysignal = mysignal * len(x1)
    
    return mysignal, Fa1
    
    
def create_artificial_ml(num_samples, mean = 1, std = 1/50):

    incr = np.arange(1600, 1700, (1700-1600)/num_samples)
    s3 = np.random.normal(mean, std, size=num_samples) * incr 
        
    return s3
        
    
