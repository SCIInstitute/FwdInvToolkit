# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 08:48:45 2016

@author: andrewcorbato
"""

def lowpassma(M,win):
    
#     function lowpassma(M,win):
#     returns LOW
#     lowpass Moving Average filter applied to rows of matrix M; finite impulse
#     response; no phase shift
#     zero output for frequencies at multiples of  f=win*sampling_interval
#     for highpass: use PSI=PSI-lowpassma(M,win)
#     2005-02-17;  %V. Jacquemet/ A. van Oosterom
#     2008-02-12; AvO winb=np.round()   changed into:   winb=np.floor() 
    
    import numpy as np
    Low = M
    if win < 2:
        return Low
    trans = 0
    dims_M = np.shape(M)
    nsig = dims_M[0]
    nt = dims_M[1]
    if nt == 1:
        trasp = 1
        M = np.transpose(M)
        dims_M = np.shape(M)
        nsig = dims_M[0]
        nt = dims_M[1]
    winb = np.floor(win/2)
    wine = win-winb
    LEAD = np.dot(M[:,0],np.ones(1,winb))
    TRAIL = np.dot(M[:,0],np.ones(1,wine))
    dims_LEAD = np.shape(LEAD)
    dims_TRAIL = np.shape(TRAIL)
    rows = dims_M[0]
    columns = dims_LEAD[1]+dims_M[1]+dims_TRAIL[1]
    storage = np.zeros((rows,columns))
    storage[:,0:dims_LEAD[1]] = LEAD
    storage[:,dims_LEAD[1]:dims_LEAD[1]+dims_M[1]] = M
    storage[:,dims_LEAD[1]+dims_M[1]:dims_LEAD[1]+dims_M[1]*dims_TRAIL[1]] = TRAIL
    M = storage
    X = np.cumsum(np.concatenate(np.zeros((nsig,1)),M),axis=1)
    LOW = X[:,win:win+nt]-X[:,0:nt]
    LOW = np.dot(LOW,np.linalg.inv(win))
    if trans == 1:
        LOW = np.transpose(LOW)
        
    return LOW
    
    
    
            