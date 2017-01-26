# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 09:45:02 2016

@author: andrewcorbato
"""

def InterpolateCurve(CurveParameters,InterpolationDensity,varargin):
    # author: Burak Erem
    # input: varargin should be a list


    import numpy as np
    
    # check type of varargin
    if type(varargin) != type([]):
        varargin = list(varargin)
        if type(varargin) != type([]):
            import sys
            sys.exit('ERROR: varargin must be a python list or compatible w/ list().')
    
    # input sizes and number of elements
    dims_CurveParameters = np.shape(CurveParameters)
    prod_dims_CurveParameters = np.prod(dims_CurveParameters,axis=0)
    dims_dims_CurveParameters = np.shape(dims_CurveParameters)
    numel_dims_dims = dims_dims_CurveParameters[0]*dims_dims_CurveParameters[1]    
    PeriodicityFlag = np.zeros((numel_dims_dims-1,1))
    numel_varg = len(varargin)
    
    if numel_varg > 0:
        p = 'periodic'
        for ii in range(0,numel_varg,1):
            if varargin[ii].lower() == p:
                # unless this is followed by a numeric array of dimension
                # indices that should be periodic,
                # assume they are all periodic
                if numel_varg > ii:
                    if type(varargin[ii+1]) != type(p):
                        PeriodicityFlag(varargin[ii+1]) = 1
                        dims_PeriodicityFlag = np.shape(PeriodicityFlag)
                        numel_PeriodicityFlag = dims_PeriodicityFlag[0]*dims_PeriodicityFlag[1]
                        if numel_PeriodicityFlag > (numel_dims - 1):
                            print('WARNING: Dimensions specified as being periodic \
                            exceed input dimensions. \n')
                else:
                    PeriodicityFlag = np.ones((dims_PeriodicityFlag[0],dims_PeriodicityFlag[1]))
    
    # Form 1-D spline interpolation matrices of appropriate sizes
    # TensorEdgeDimensions = np.sort(np.unique(dims_CurveParameters[2,-1]))
    # for jj = TensorEdgeDimensions
    #   SplineMatrix[i-2]=np.transpose(splineMatrix(i-2,1,InterpolationDensity))

    
    TensorEdgeDimensions = dims_CurveParameters[1:-1]
    dims_TensorEdgeDimensions = np.shape(TensorEdgeDimensions)
    numel_Tensor = dims_TensorEdgeDimensions[0]*dims_TensorEdgeDimensions[1]
    SplineMatrix = [None]*numel_Tensor
    for jj in range(0,numel_Tensor,1):
        if PeriodicityFlag[jj] == 0: # if not periodic
            SplineMatrix[jj] = np.transpose(splinematrix(TensorEdgeDimensions[jj]-2,1,InterpolationDensity))
        else: # if periodic
            SplineMatrix[jj] = np.transpose(splinematrix(TensorEdgeDimensions[jj]-2,1,InterpolationDensity,'Periodic'))
    
    # intialize interpolated curves as curve parameters
    InterpolatedCurve = CurveParameters
    
    # Interpolate spline curves using tensor-matrix "right" multiplication, for
    # all the "right-hand" sides (i.e. tensor indices, not including the first one)
    for kk in range(0,(numel_dims_dims-1),1):
        InterpolatedCurve = tensorRightMatrixMultiply(InterpolatedCurve,kk,SplineMatrix[kk])
    
    return InterpolatedCurve
        
    