# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 09:34:55 2016

@author: andrewcorbato
"""

def initializeCurveParamsFromTimeSeries(Data,NumberOfKnots):
    # author: Burak Eram
    import numpy as np
    dims_Data = np.shape(Data)
    T = dims_Data[1]
    
    knotinds = np.floor(np.linspace(0,T,NumberOfKnots))
    Knots = Data[:,knotinds]
    
    FirstKnotDeriv = (Knots[:,1]-Knots[:,0])/(knotinds(1)-knotinds(0))
    LastKnotDeriv = (Knots[:,-1]-Knots[:,-2])/(knotinds(-1)-knotinds(-2))
    
    dims_FirstKnotDeriv = np.shape(FirstKnotDeriv)
    dims_LastKnotDeriv = np.shape(LastKnotDeriv)
    dims_Knots = np.shape(Knots)
    columns = dims_FirstKnotDeriv[1]+dims_Knots[1]+dims_LastKnotDeriv[1]
    CurveParams = np.zeros((dims_Knots[0],columns))
    CurveParams[:,0:dims_FirstKnotDeriv[1]] = FirstKnotDeriv
    CurveParams[:,dims_FirstKnotDeriv[1]:dims_FirstKnotDeriv[1]+dims_Knots[1]] = Knots
    CurveParams[:,dims_FirstKnotDeriv[1]+dims_Knots[1]:columns] = LastKnotDeriv
    
    return CurveParams