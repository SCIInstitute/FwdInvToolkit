# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:22:17 2016

@author: andrewcorbato
"""
def icp(model,data,*arg):
#   ICP Iterative Closest Point Algorithm. Takes use of
#   Delaunay tesselation of points in model.
#
#   Ordinary usage:
#   [R, S, T,out] = icp(model,data)
#
#   ICP fit points in data to the points in model.
#   errors with the closest model points and data points.
#
#   INPUT:
#   model - matrix with model points, [Pm_1 Pm_2 ... Pm_nmod]
#   data - matrix with data points,   [Pd_1 Pd_2 ... Pd_ndat]
#
#   OUTPUT:
#   R - rotation matrix and
#   S - scaling value (uniform)
#   T - translation vector accordingly so
#
#           newdata = R*data*S + T .
#   out - output parameters of the algorithm
#       out.iter - number of iterations used
#       out.error - final error of the algorithm
#
#   newdata are transformed data points to fit model
#
#
#   Special usage:
#   icp(model)  or icp(model,tes_flag)
#
#   ICP creates a Delaunay tessellation of points in
#   model and save it as global variable Tes. ICP also
#   saves two global variables ir and jc for tes_flag=1 (default) or
#   Tesind and Tesver for tes_flag=2, which
#   makes it easy to find in the tesselation. To use the global variables
#   in icp, put tes_flag to 0.
#
#   Other usage:
#   [R, S, T] = icp(model,data,max_iter,min_iter,...
#                         fitting,thres,init_flag,tes_flag,delta)
#
#   options: if you would like to change options for the default call in\
# function inputs as dictionary, copy and paste the following and edit if desired:
# 
#       options = {'max_iter': 104, 
#       'min_iter': 4,
#       'fitting': 2,
#       'thres': 1e-05,
#       'init_flag': 1,
#       'tes_flag': 1,
#       'delta': 10,
#       'mode': 'rigid',
#       'refpnt': np.array([]),}
#
# 	max_iter - maximum number of iterations. Default=104
#
# 	min_iter - minimum number of iterations. Default=4
#
#   fitting  -  =2 Fit with respect to minimize the sum of square errors. (default)
#                  alt. =[2,w], where w is a weight vector corresponding to data.
#                  w is a vector of same length as data.
#                  Fit with respect to minimize the weighted sum of square errors.
#               =3 Fit with respect to minimize the sum to the amount 0.95
#                  of the closest square errors.
#                  alt. =[3,lambda], 0.0<lambda<=1.0, (lambda=0.95 default)
#                  In each iteration only the amount lambda of the closest
#                  points will affect the translation and rotation.
#                  If 1<lambda<=size(data,2), lambda integer, only the number lambda
#                  of the closest points will affect the translation and
#                  rotation in each iteration.
#
# 	thres - error differens threshold for stop iterations. Default 1e-5
#
# 	init_flag  -  =0 no initial starting transformation
#                 =1 transform data so the mean value of
#                     data is equal to mean value of model.
#                     No rotation. (init_flag=1 default)
#
# 	tes_flag  -  =0 No new tesselation has to be done. There
#                   alredy exists one for the current model points.
#                =1 A new tesselation of the model points will
#                   be done. (default)
#                =2 A new tesselation of the model points will
#                   be done. Another search strategy than tes_flag=1
#                =3 The closest point will be find by testing
#                   all combinations. No Delaunay tesselation will be done.
#
#   delta -      > 0  scaling parameter sensitivity.  this option allows
#                     restriction of the scaling parameter to prevent the
#                     algorithm from converging to a single point.  delta
#                     of 0 produces no scaling changes, 1 can double the
#                     scale. 'rigid' mode only.
#
#   mode - registration method used.  There are currently two options:
#                 'affine' and 'rigid'.  default is 'rigid'.
#
#
#   refpnt - (optional) (An empty vector is default.) refpnt is a point corresponding to the
#                 set of model points wich correspondig data point has to be find.
#                 How the points are weighted depends on the output from the
#                 function weightfcn found in the end of this m-file. The input in weightfcn is the
#                 distance between the closest model point and refpnt.
#
#   To clear old global tesselation variables run: "clear global Tes ir jc" (tes_flag=1)
#   or run: "clear global Tes Tesind Tesver" (tes_flag=2) in Command Window.
#
#   based on an m-file can be downloaded for free at
#   http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=12627&objectType=FILE
#
#

    import numpy as np
    import scipy as sp
    from numpy import matlib as ml
    from numpy import matrix as mx
    import sys

    if model.any() == False and data.any() == False:
        raise ValueError('ERROR: Only one input. There must be at least two inputs.')
    elif model.any() == False:
        raise ValueError('ERROR: Model cannot be an empty matrix.')
    global MODEL, DATA, TR, TT, TS
    
    dims_model = np.shape(model)
    if dims_model[1]<dims_model[0]:
        MODEL = np.transpose(model)
    else:
        MODEL = model
        
    dims_data = np.shape(data)
    if dims_data[1]<dims_data[0]:
        data = np.transpose(data)
        dims_data = np.shape(data)
        DATA = data
    else:
        DATA = data
        
    if dims_data[0] is not dims_model[0]:
        raise ValueError('ERROR: DATA and MODEL cannot have different dimensions.')
# default options
    if not len(arg)==0:
        opt = arg[0]
    else:
        opt = {'max_iter': 104, 'min_iter': 4,'fitting': [2],'thres': 1e-05, \
        'init_flag': 1,'tes_flag': 1, 'delta': 10,'mode': 'rigid', \
        'refpnt': np.array([]),}

    print(opt)

    if not isinstance(opt, dict):
        raise ValueError("opt must be dictionary of options")
    if not 'max_iter' in opt:
        opt['max_iter'] = 104;
    if not 'min_iter' in opt:
        opt['min_iter'] = 4;
    if not 'fitting' in opt:
        opt['fitting'] = [2];
    if not 'thres' in opt:
        opt['thres'] = 1e-5;
    if not 'delta' in opt:
        opt['delta'] = 10;
    if not 'mode' in opt:
        opt['mode'] = 'rigid';
    if not 'refpnt' in opt:
        opt['refpnt'] = np.array([]);
    if not 'init_flag' in opt:
        opt['init_flag'] = 1;
    if not 'tes_flag' in opt:
        opt['tes_flag'] = 1;




# move options out of dictionary for ease of use
    max_iter = opt['max_iter']
    min_iter = opt['min_iter']
    fitting = opt['fitting']
    thres = opt['thres']
    init_flag = opt['init_flag']
    tes_flag = opt['tes_flag']
    delta = opt['delta']
    mode = opt['mode']
    refpnt = opt['refpnt']
    
# input error checks
    if (tes_flag != 0 and tes_flag != 1 and tes_flag != 2 and tes_flag != 3):
        tes_flag = 3
        print('WARNING: tes_flag has been set to the default value of 3, \
must be 0-3.')    
    if dims_model[1] == dims_model[0] and tes_flag is not 0:
        import sys
        raise ValueError('ERROR: This icp method requires the number of model points \
to be greater than the dimension')
    if max_iter < min_iter:
        max_iter = min_iter
        print('WARNING: max_iter has been set to equal min_iter.')
    if min_iter < 0:
        min_iter = 0
        print('WARNING: min_iter has been set to zero.')
    if thres < 0:
        thres= np.abs(thres)
        print('WARNING: thres negative, has been set to its absolute value.')
    if fitting[0] == 2:
        dims_fitting = np.shape(fitting)
        lef = np.max(dims_fitting)
        if lef > 1:
            if dims_fitting[0] < dims_fitting[1]:
                fitting = np.transpose(fitting)
            if lef < (dims_data[1]+1):
                print('WARNING: Illegal size of fitting; unweighted minimization will be used.')
                fitting = 2
            elif np.min(fitting[1:dims_data[1]+1]) < 0:
                print('WARNING: Illegal value of the weights; unweighted minimization will be used.')
                fitting = 2
            elif np.max(fitting[1:dims_data[1]+1]) == 0:
                fitting = 2
                print('WARNING: Illegal values of the weights; unweighted minimization will be used.')
            else:
                su = np.sum(fitting[1:dims_data[1]+1])
                fitting[1:dims_data[1]+1] = fitting[1:dims_data[1]+1]/su
    elif fitting[0] == 3:
        if np.max(dims_fitting) < 2:
            fitting = np.array(fitting,np.rint(0.95*dims_data[1]))
        elif fitting[1] > 1:
            if fitting[1] > np.floor(fitting[1]):
                fitting[1] = np.floor(fitting[1])
                f = fitting[1] 
                print('WARNING: lambda has been set to: %i' %f)
            if fitting[1] > dims_data[1]:
                fitting[1] = dims_data[1]
                f = fitting[1]
                print('WARNING: lambda has been set to: %i' %f)
        elif fitting[1] > 0:
            if fitting[1] <= 0.5:
                print('WARNING: lambda is small problems may arise')
            fitting[1] = np.rint(fitting[2]*dims_data[1])
        elif fitting[1] <= 0:
            fitting[1] = np.rint(0.95*dims_data[1])
            f = fitting[1]
            print('WARNING: lambda has been set to: %i' %f)
    else:
        fitting = 2
        print('WARNING: fitting has been reset to its default value of 2')
    if init_flag is not 0 or init_flag is not 1:
        init_flag = 1
        print('WARNING: init_flag has been reset to its default value of 1')
    dims_refpnt = np.shape(refpnt)
    if dims_refpnt[1] > dims_refpnt[0]:
        refpnt = np.transpose(refpnt)
    if dims_refpnt[0] is not dims_data[0]:
        refpnt = np.array([])
        print('WARNING: Dimensions of refpnt do not match data, refpnt has \
been reset to its default value of an empty array')
    if dims_refpnt[1]>1:
        refpnt = refpnt[:,1]
        print('WARNING: Only the first point in refpnt has been used.')
        
#   start the ICP algorithm
    N = dims_data[1]
    icp_init(init_flag,fitting) # initiate a starting rotation matrix\
# and starting translation vector
    tes_flag = icp_struct(tes_flag) # construt a Delauney tesselation and\
# two vectors to make it easy to find in Tes
    ERROR = icp_closest_start(tes_flag,fitting) # initiate a vector with\
# indices of closest MODEL points
    icp_transformation(fitting,np.array([]),delta,mode) # find transformation
    DATA = TR*data*TS
    DATA = DATA+ml.repmat(TT,1,N)
    for iter in range(0,max_iter,1):
        olderror = ERROR
        ERROR = icp_closeste(tes_flag,fitting)
        if iter < min_iter:
            icp_transformation(fitting,np.array([]),delta,mode)
        else:
            icp_transformation(fitting,refpnt,delta,mode)
        
        DATA = TR*data*TS
        DATA = DATA+ml.repmat(TT,1,N)
        if iter >= min_iter:
            if np.abs(olderror-ERROR) < thres:
                break
            
    class out():
        def __init__(self):
            self.iter = iter
        def __init__(self):
            self.error = ERROR
        def __init__(self):
            self.residual = olderror-ERROR
    out = out()
    
    return TR, TS, TT, out
###############################################################################   
def icp_init(init_flag,fitting):
    # function icp_init(init_flag)
    # ICP_INIT intial alignment for ICP.
    global MODEL, DATA, TR, TT, TS
    if init_flag == 0:
        dims_MODEL = np.shape(MODEL)        
        TR = np.eye(dims_MODEL[0],dtype=int)
        TT = np.zeros(dims_model[0],0)
        TS = 1
    elif init_flag ==1:
        dims_DATA = np.shape(DATA)
        N = dims_DATA[1]
        if fitting[0] == 2:
            dims_fitting = np.shape(fitting)
            if mx.size(fitting) == 1:
                mem = np.mean(MODEL,1)
                med = np.mean(DATA,1)
            else:
                mem = np.mean(MODEL,1)
                med = DATA*fitting[1:(N)]
        else:
            mem = np.mean(MODEL,1)
            med = np.mean(DATA,1)
        TR = np.eye(dims_MODEL[0],dtype=int)
        TT = mem-med
        TS = 1
        DATA = DATA + ml.repmat(TT,1,N) # apply transformation
    else:
        import sys
        raise ValueError('ERROR: Wrong init_flag')
###############################################################################    
def icp_struct(tes_flag):
    global Tes, ir, jc, Tesind, Tesver
    if tes_flag != 3:
        if tes_flag == 0:
            global ir
            dims_ir = np.shape(ir)
            if dims_ir[0] == 0:
                global Tesind
                dims_Tesind = np.shape(Tesind)
                if dims_Tesind[0] == 0:
                    import sys
                    raise ValueError('ERROR: No tesselation system exists')
                else:
                    tes_flag = 2
            else:
                tes_flag = 1
    elif tes_flag == 3:
        return tes_flag
    else:
      
        [m,n] = np.shape(MODEL)
        if m == 1:
            ind1 = np.argsort(MODEL)
            so1 = np.sort(MODEL) 
            Tes = np.zeros((n-1,1))
            Tes[0:(n-2),0] = np.transpose(ind1[1:n])
            Tes[1:n-1,1] =np.transpose(ind1[1:(n-1)])               
            Tes[0,1] = Tes[0,0]
            Tes[n-1,0] = Tes[n-1,1]
            Tes[ind1,:] = Tes[:,:]
        else:
            Tes = delaunayn(np.transpose(MODEL))
            dims_Tes = np.shape(Tes)
            if dims_Tes[0] == 0:
                mem = np.mean(MODEL,axis=1)
                MODELm = MODEL-ml.repmat(mem,0,n)
                [U, S, V] = svd(MODELm*np.transpose(MODELm))
                onbasT = np.transpose(U[:,(np.diag(S)>1e-08)])
                MODELred = onbasT*MODEL
                dims_MODELred = np.shape(MODELred)
                if dims_MODELred[0] == 1:
                    ind1 = np.argsort(MODELred)
                    so1 = np.sort(MODELred)
                    Tes = np.zeros((n,2))
                    Tes[0:(n-2),0] = np.transpose(ind1[1:n-1])
                    Tes[1:n-1,1] =np.transpose(ind1[0:(n-2)])               
                    Tes[0,1] = Tes[0,0]
                    Tes[n-1,0] = Tes[n-1,1]
                    Tes[ind1,:] = Tes[:,:]
                else:
                    Tes = delaunyan(np.transpose(MODELred))
    Tes = np.sort(np.sort(Tes,axis=1),axis=1)
    [mT,nT] = np.shape(Tes)
    if tes_flag == 1:
        num = np.zeros((1,n))
        for jj in range(0,mT,1):
            for kk in range(0,nT,1):
                num[Tes[jj,kk]] = num[Tes[jj,kk]]+1
        num = np.cumsum(num,axis=0)
        jc = np.ones(1,n+1)
        jc[1:-1] = num + jc[1:-1]
        ir = np.zeros((1,num[-1]))
        ind = jc[0:-2]
        
        # calculate ir
        for i in range(0,mT,1):
            for j in range(0,nT,1):
                ir[ind[Tes[i,j]]] = i
                ind[Tes[i,j]] = ind[Tes[i,j]]+1
    else: # tes_flag ==2
        Tesind = np.zeros(mT,nT)
        Tesver = np.zeros(mT,nT)
        couvec = np.zeros(mT,1)
        
        for i in range(0,mT-1,1):
            for j in range(i+1,mT,1):
                if couvec[i] == nT:
                    break
                elif Tes[i,0] == Tes[j,0]:
                    if nT == 2:
                        Tesind[i,1] = j
                        Tesind[j,1] = i
                        Tesver[i,1] = 2
                        Tesver[j,1] = 2
                        couvec[i] += 1
                        couvec[j] += 1
                    else:
                        for k in range(1,nT,1):
                            for kk in range(k,nT,1):
                                if np.all(np.concatenate(Tes[i,1:k-1],Tes[i,k+1:nT-1]) == np.concatenate(Tes[j,1:k-1],Tes[j,k+1:nT-1])):
                                    Tesind[i,k] = j
                                    Tesind[j,kk] = i
                                    Tesver[i,k] = kk
                                    Tesver[j,kk] = k
                                    couvec[i] +=1
                                    couvec[j] +=1   
                            if couvec[i] == nT or couvec[j] == nT:
                                break
                elif Tes[i,0] == Tes[j,1]:
                    if nT==2:
                        Tesind[i,1] = j
                        Tesind[j,0] = i
                        Tesver[i,1] = 1
                        Tesver[j,0] = 2
                        couvec[i] +=1
                        couvec[j] +=1
                    else:
                        for k in range(1,nT,1):
                            if np.all(np.concatenate(Tes[i,1:k-1],Tes[i,k+1:nT-1]),axis=0) == Tes[j,2:nT-1]:
                                Tesind[i,k] = j
                                Tesind[j,k] = i
                                Tesver[i,k] = 1
                                Tesver[j,0] = k
                                couvec[i] +=1
                                couvec[j] +=1
                            if couvec[i] == nT or couvec[j] == nT:
                                break
                elif Tes[i,1] == Tes[j,0]:
                    if nT == 2:
                        Tesind[i,0] = j
                        Tesind[j,0] = i
                        Tesver[ui0] = 1
                        Tesver[j,0] = 1
                        couvec[i] +=1
                        couvec[j] +=1
                        if Tes[j,0] > Tes[i,1]:
                            break
                    elif np.all(Tes[yy,2:-1]) == Tes[zz,2:-1]:
                        Tesind[i,0] = j
                        Tesind[j,0] = i
                        Tesver[i,0] = 1
                        Tesver[j,0] = 1
                        couvec[i] +=1
                        couvec[j] +=1
                elif Tes[j,0] > Tes[i,1]:
                    break
    return tes_flag
###############################################################################    
def icp_closest_start(tes_flag,fitting):
# Usage:
#       ERROR = icp_closest(tes_flag)
# ERROR=sum of all errors between points in MODEL and points in DATA.
#
# ICP_CLOSEST_START finds indexes of closest MODEL points for each point in DATA.
# The _start version allocates memory for iclosest and finds the closest MODEL points
# to the DATA points
    global MODEL, DATA, Tes, ir, jc, iclosest
    if tes_flag == 3:
        dims_MODEL = np.shape(MODEL)
        dims_DATA = np.shape(DATA)
        mm = dims_MODEL[1]
        md = dims_DATA[1]
        
        iclosest = np.zeros((1,md))
        ERROR = 0;
        for ID in range(0,md,1):
            dist = float('inf')
            for im in range(0,mm,1):
                dista = np.linalg.norm((MODEL[:,im]-DATA[:,ID]))
                if dista < dist:
                    iclosest[ID]=im
                    dist = dista
        ERROR += err(dist,fitting,ID)
        
    elif tes_flag == 1:
        dims_DATA = np.shape(DATA)
        md = dims_DATA[1]
        dims_MODEL = np.shape(MODEL)
        iclosest = np.zeros((1,md))
        mid = np.round(md/2)
        iclosest[mid] = np.round((dims_MODEL[1]/2))
        bol = 1
        while bol:
            bol = np.logical_not(bol)
            distc = np.linalg.norm((DATA[:,mid]-MODEL[:,iclosest[mid]]))
            distcc = 2*distc
            for i in range(ir[jc[iclosest[mid]]],ir[jc[iclosest[mid]+1]-1],1):
                for ind in Tes[i,:]:
                    distcc = np.linalg.norm((DATA[:,mid]-MODEL[:,ind]))
                    if distcc<distc:
                        distc = distcc
                        bol = np.logical_not(bol)
                        iclosest[mid] = ind
                        break
                if bol:
                    break
        ERROR = err(distc,fitting,mid)
    
        for ID in range((mid+1),md,1):
            iclosestID = iclosest[ID-1]
            bol = np.logical_not(bol)
            while bol:
                bol = np.logical_not(bol)
                distc = np.linalg.norm((DATA[:,ID]-MODEL[:,iclosest[ID]]))
                distcc = 2*distc
                for i in range(ir[jc[iclosest[ID]]],ir[jc[iclosest[ID]+1]-1],1):
                    for ind in Tes[i,:]:
                        distcc = np.linalg.norm((DATA[:,ID]-MODEL[:,ind]))
                        if distcc<distc:
                            distc = distcc
                            bol = np.logical_not(bol)
                            iclosest[mid] = ind
                            break
                    if bol:
                        break
            ERROR += err(distc,fitting,ID)
    
        for ID in range((mid-1),-1,1):
            iclosest[ID] = iclosest[ID+1]
            bol = np.logical_not(bol)
            while bol:
                bol = np.logical_not(bol)
                distc = np.linalg.norm((DATA[:,ID]-MODEL[:,iclosest[ID]]))
                distcc = 2*distc
                for i in range(ir[jc[iclosest[ID]]],ir[jc[iclosest[ID]+1]-1],1):
                    for ind in Tes[i,:]:
                        distcc = np.linalg.norm((DATA[:,ID]-MODEL[:,ind]))
                        if distcc<distc:
                            distc = distcc
                            bol = np.logical_not(bol)
                            iclosest[mid] = ind
                            break
    else: # tes_flag == 2
        dims_DATA = np.shape(DATA)
        md = dims_DATA[1]
        iclosest = np.zeros((1,md))
        icTesind = np.zeros((1,md))
        dims_Tes = np.shape(Tes)
        mTes = dims_Tes[0]
        nTes = dims_Tes[1]
        mid = np.round(md/2)
        icTesind[mid] = np.round(mTes/2)
        iclosest[mid] = np.max(Tesind[icTesind[mid],:])
        visited = np.zeros(1,mTes)
        visited[icTesind[mid]] = 1
        presum = (ml.repmat(DATA[:,mid],1,nTes)-MODEL[:,Tes[icTesind[mid],:]])**2
        postsum = np.sum(presum,0)
        di2vec = np.sqrt(postsum)
        bol = 1
        while bol:
            ind = np.argsort(di2vec)
            so = np.sort(di2vec)
            for ii in range(nTes,-1,2):
                Ti = Tesind[icTesind[mid],ind[ii]]
                if Ti > 0:
                    if np.logical_not(visited(Ti)):
                        break
            if Ti == 0:
                bol = np.logical_not(bol)
            elif visited(Ti):
                bol = np.logical_not(bol)
            else:
                Tv = Tesver[icTesind[mid],ind[ii]]
                visited[Ti] = 1
                icTesind[mid] = Ti
                di2vec[1:(Tv-1),(Tv+1):nTes] = di2vec[1:(ind[ii]-1),(ind[ii]+1):nTes]
                di2vec[Tv] = np.linalg.norm((DATA[:,mid]-MODEL[:,Tes(Ti,Tv)]))
        iclosest[mid] = Tes[icTesind[mid],ind[1]]
        ERROR = err(so[1],fitting,mid)
        
        for ID in range((mid+1),md,1):
            iclosest[ig] = iclosest[ID-1]
            icTesind[ig] = icTesind[ID-1]
            visited[:] = 0
            visited[icTesind[ID-1]] = 1
            presum = (ml.repmat(DATA[:,ID],1,nTes)-MODEL[:,Tes[icTesind[ID],:]])**2
            postsum = np.sum(presum,0)
            di2vec = np.sqrt(postsum)
            bol = np.logical_not(bol)
            while bol:
                so = np.sort(di2vec)
                ind = np.argsort(di2vec)
                for ii in range(nTes,-1,2):
                    Ti = Tesind[icTesind[ID]],ind[ii]
                    if Ti > 0:
                        if np.logical_not(visited[Ti]):
                            break
                if Ti == 0:
                    bol = np.logical_not(bol)
                elif visited[Ti]:
                    bol = np.logical_not(bol)
                else:
                    Tv = Tesver[icTesind[ID],ind[ii]]
                    visited[Ti] = 1
                    icTesind[ig] = Ti
                    di2vec[1:(Tv-1),(Tv+1):nTes] = di2vec[1:(ind[ii]-1),(ind[ii]+1):nTes]
                    di2vec[Tv] = np.linalg.norm((DATA[:,mid]-MODEL[:,Tes(Ti,Tv)]))
            iclosest[ID] = Tes[icTesind[ID],ind[1]]
            ERROR += err(so[1],fiiting,ID)
    return ERROR
###############################################################################               
def icp_transformation(fitting,refpnt,delta,optmode):
    
    global MODEL, DATA, iclosest, TR, TT, TS
    dims_DATA = np.shape(DATA)
    M = dims_data[0]
    N = dims_data[1]
    bol = 0
    if options['refpnt'] == True:
        bol = 1
        dis = np.sqrt(np.sum((MODEL[:,iclosest]-ml.repmat(refpnt,1,N))**2))
        weights = weightfcn(np.transpose(dis))
    if bol:
        if fitting[0] == 2:
            dims_fitting = np.shape(fitting)
            if dims_fitting[0] > 1 or dims_fitting[1] > 1:
                weights = fitting[1:(N+1)]*weights
                weights = weights/np.sum(weights,axis=0)
            med = DATA*weights
            mem = MODEL[:,iclosest]*weights
            A = DATA - ml.repmat(med,1,N)
            B = MODEL[:,iclosest] - ml.repmat(meme,1,N)
            
            for ii in range(0,N,1):
                A[:,ii] = weights[ii]*A[:,ii]
        elif fitting[1] == 3:
            V = np.sum((DATA-MODEL[:,iclosest])**2,axis = 0)
            ind = np.argsort(V)
            soV = np.sort(V)
            ind = ind[0:fitting[1]-1]
            weights[ind] = weights[ind]/np.sum(weights[ind])
            med = DATA[:,ind]*weihgts[ind]
            mem = MODEL[:,iclosest[ind]]*weights[ind]
            A = DATA[:,ind]-ml.repmat(mem,1,fitting[1])
            B = MODEL[:,iclosest[ind]]-ml.repmat(mem,1,fitting[1])
            
            for ii in range(0,fitting[1],1):
                A[:,ii] = weights[ind[ii]]*A[:,ind[ii]]
    else:
        if fitting[0] == 2:
            dims_fitting = np.shape(fitting)
            if dims_fitting[0] > 1 or dims_fitting[1] > 1:
                if mode == 'affine':
                    P = np.zeros((3*N,12))
                    q = np.zeros((3*N,12))
                    ind = np.array([1,2,3])
                    for n in range(0,N,1):
                        a = np.concatenate(DATA[0,n],0,0,DATA[1,n],0,0,DATA[2,n],0,0,1,0,0, axis=1)
                        b = np.concatenate(0,DATA[0,n],0,0,DATA[1,n],0,0,DATA[2,n],0,0,1,0, axis=1)
                        c = np.concatenate(0,0,DATA[0,n],0,0,DATA[1,n],0,0,DATA[2,n],0,0,1, axis=1)
                        P[ind,:] = np.concatenate(a,b,c, axis=0)
                        q[ind] = MODEL[0:2,iclosest[n]]
                        ind = ind+3
                elif mode == 'rigid':
                    med = np.mean(DATA,axis=1)
                    mem = MODEL[:,iclosest]*fitting[1:(N+1)]
                    A = DATA - ml.repmat(med,1,N)
                    B = MODEL[:,iclosest] - ml.repmat(mem,1,N)
            else:
                med = DATA*fitting[1:(N+1)]
                mem = MODEL[:,iclosest]*fitting[1:(N+1)]
                A = DATA - ml.repmat(med,1,N)
                B = MODEL[:,iclosest] - ml.repmat(mem,1,N)
                for ii in range(0,N,1):
                    A[:,ii] = fitting[ii+1]*A[:,ii]
        elif fitting[1] == 3:
            V = np.sum((DATA-MODEL[:,iclosest])**2,axis=0)
            ind = np.argsort(V)
            soV = np.sort(V)
            ind = ind[0:fitting[1]-1]
            
            if mode == 'rigid':
                med = np.mean(DATA[:,ind],axis=1)
                mem = np.mean(MODEL[:,iclosest[ind]],axis=1)
                A = DATA[:,ind] - ml.repmat(med,1,fitting[1])
                B = MODEL[:,iclosest[ind]] - mlrepmat(mem,1,fitting[1])
                
            elif mode == 'affine':
                p_ind = ind
                P = np.zeros((3*fitting[1],12))
                q = np.zeros((3*fitting[1],1))
                ind = np.array([1,2,3])
                for n in range(0,fitting[1],1):
                    a = np.concatenate(DATA[0,n],0,0,DATA[1,n],0,0,DATA[2,n],0,0,1,0,0, axis=1)
                    b = np.concatenate(0,DATA[0,n],0,0,DATA[1,n],0,0,DATA[2,n],0,0,1,0, axis=1)
                    c = np.concatenate(0,0,DATA[0,n],0,0,DATA[1,n],0,0,DATA[2,n],0,0,1, axis=1)
                    P[ind,:] = np.concatenate(a,b,c, axis=0)
                    q[ind] = MODEL[0:2,iclosest[p_ind[n]]]
                    ind = 3
    if mode == 'affine':
        theta = q/P
        a_r = np.concatenate(theta[0],theta[3],theta[6],axis=0)
        b_r = np.concatenate(theta[1],theta[4],theta[7],axis=0)
        c_r = np.concatenate(theta[2],theta[5],theta[8],axis=0)
        R = np.concatenate(a_r,b_r,c_r,axis=1)
        T = np.concatenate(theta[9],theta[10],theta[11],axis=1)
        S = 1
    elif mode == 'rigid':
        normA = np.sqrt(np.trace(A*np.transpose(A)))
        normB = np.sqrt(np.trace(B*np.transpose(B)))
        A = A/normA
        B = B/normB
        U, Sig, V = np.linalg.svd(np.multiply(B,np.transpose(A)))
        U[:,-1] = U[:,-1]*np.linalg.determinant(U*np.transpose(V))
        R = U*np.transpose(V)
        trsqrtAB = np.sum(np.diag(Sig))
        S = trsqrtAB*normB/normA
        max_scale = 1+delta
        min_scale = 1/max_scale
        if (TS*S)>max_scale:
            S = max_scale/TS
        if (TS*S)<min_scale:
            S = min_scale/TS
        # compute the translation
        T = (mem-R*S*med)
    TR = R*TR
    TT = R*TT*S+T
    TS = S*TS
###############################################################################        
def err(dist,fitting,ind):
    if fitting[1] == 2:
        dims_fitting = np.shape(fitting)
        if (ind+1) > dism_fitting[0] and (ind+1) > dism_fitting[1]:
            ERR = dist**2
        else:
            ERR = fitting[ind+1]*dist**2
    elif fitting[1] == 3:
        ERR = dist**2
    else:
        ERR = 0
        print('WARNING: Unknown fitting value.')
    return ERR
###############################################################################        
def weightfcn(distances):
    max_distances = np.amax(distances)
    min_distances = np.amin(distances)
    if max_distances > (1.1*min_distances):
        weights=1+min_distances/(max_distances-min_distances)-distances/(max_distances-min_distances)
    else:
        weights = max_distances+min_distances-distances
    weights = weights/np.sum(weights)
    return weights
