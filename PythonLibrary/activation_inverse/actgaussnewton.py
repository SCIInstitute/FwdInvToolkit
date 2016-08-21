# -*- coding: utf-8 -*-
def ActGaussNewton(A, Y, L, tauinit, Lambda, w, minstep):
#   Implements the Gauss-Newton algorithm for solving the activation-based
#   inverse problem of electrocardiography.
#   => minimizes the objective function ||Y-A*X||^2+lambda*||L*X||^2 where
#   X is parameterized by the C^1 polynomial approximation to a step function
#   as explained in "The Depolarization Sequence of the Human Heart Surface
#   Computed from Measured Body Surface Potentials" by Geertjan Huiskamp and
#   Adriaan van Oosterom.
# 
#   Input Variables:
#   A: Forward matrix
#   Y: Observations (columns index time from 1 to T=size(Y,2))
#   L: Regularization matrix (typically a surface Laplacian approximation)
#   Lambda: Regularization parameter
#   w: Width parameter in step function approximation
#   tauinit: Initial phase shifts for starting the algorithm
# 
#   Output Variables:
#   tau: Solution phase shifts of the step functions

#   import scipy and numpy packages
    import scipy as sp
    import numpy as np
    
    L = sp.sparse.csc_matrix.todense(L)
    tau = tauinit()
    dims_tau = np.shape(tauinit)
    step = 2*np.ones(dims_tau)
#   num = linesteps the number of steps in [0,1] to consider for the line search
    alpha = np.linspace(0,1, num = 100)
    dims_alpha = np.shape(alpha)
    dims_N = np.shape(A)
    N = dims_N[1]
    dims_T = np.shape(Y)
    T = dims_T[1]
    u = np.linspace(1,T,T)

#   if width variable isn't an array, make it a constant array
# need reshape    dims_w = np.shape(w)
    if dims_w[1] == 1:
        w = w * np.ones((N,1)) 
        Iter = 0
        
    while np.linalg.norm(step)> minstep:
#       calculate the jacobian matrix and the residuals
        J = agnjacobian(A,Y,L,tau,Lambda,w)
        r = agnresidual(A,Y,L,tau,Lambda,w)
    
#       calculate the step direction
        G = np.multiply(np.transpose(J),J)
        size_G = np.shape(G)
        while np.cond(G)<=(np.info(float).eps):
        # possibly reciprocal of np.linalg.cond or just np.linalg.cond
#       see actgaussnewton.m line 34
            G = G+Lambda*np.identity(size_G)

        step = -G*np.linalg.inv(np.transpose(J))*r

#       perform a line search in the step direction
        err = np.zeros(dims_alpha)
        for i in range(0,dims_alpha[1],1):
            H = np.zeros((N,T))
            for n in range(1,N,1):
                H[n-1,:] = polyactrow(u-alpha[i]*step[n]-tau[n],w[n])
                
            err[i] = np.linalg.norm(Y-A*H,'fro')+ \
            Lambda*np.linalg.norm(L*H,'fro')^2
            if i>1:
                if err[i]>err[i-1]:
                    err = err[0:i]
                    break;
           
#       update step with the result of the line search
        dims_err = np.shape(err)
        vec_err = np.transpose(np.reshape(err,(dims_err[0]*dims_err[1])))
        ind = np.argmin(vec_err)
        step = alpha[ind]*step
        
#       update tau with the result of the step
        tau = tau+step
        
#       display the progress of the optimization routine
        Iter = Iter+1
        print('Step: %i\nSize: %f'%(Iter,np.linalg.norm(step)))
        return tau
###############################################################################        
def agnresidual(A,Y,L,tau,Lambda,w):

    dim_A = np.shape(A)
    dim_Y = np.shape(Y)
    N = dim_A[1]
    T = dim_Y[1]
    H = np.zeros((N,T)
    for i in range(0,N,N):
        H[i,:] = polyactrow(u-tau[i],w(i))
    E = Y-A*H
    R = np.sqrt(Lambda)*L*H
    dims_E = np.shape(E)
    dims_R = np.shape(R)
    vec_E = np.transpose(np.reshape(E,(dims_E[0]*dims_E[1])))
    vec_R = np.transpose(np.reshape(R,(dims_R[0]*dims_R[1])))
    r = np.concatenate((vec_E,vec_R), axis=0)
    return r
###############################################################################   
def agnjacobian(A,Y,L,tau,Lambda,w):
    from numpy import matlib as ml        
    dim_A = np.shape(A)
    dim_Y = np.shape(Y)
    M = dim_A[0]
    N = dim_A[1]
    T = dim_Y[1]
    dE = np.zeros((M,T,N))
    dR = np.zeros((N,T,N))
    HdH = np.zeros((1,T,N))
    u = np.linspace(1,T,T)
    for i in range(0,N,N):
        HdH[1,:,i] = dpolyactrow(u-tau[i],w(i))
    A = np.reshape(A,(M,1,N)) # MxTxL
    L = np.reshape(L,(N,1,N)) # NxTxL
    HdHm = np.ml.repmat(HdH,[M,1,1]) # MxTxL
    HdHn = np.ml.repmat(HdH,[N,1,1]) # NxTxL
    dE = np.multiply(A,HdHm)
    dR = -np.sqrt(Lambda)*np.multiply(L,HdHn)
    J = zeros(((M+N)*T,N))
    for i in range(1,N,N):
        tempE = dE[:,:,i]
        tempR = dR[:,:,i]
        dims_tempE = np.shape(tempE)
        dims_tempR = np.shape(tempR)
        vec_tempE = np.transpose(np.reshape(tempE,(dims_tempE[0]*dims_tempE[1])))
        vec_tempR = np.transpose(np.reshape(tempR,(dims_tempR[0]*dims_tempR[1])))
        J[:,i] = np.concatenate((vec_tempE,vec_tempR), axis=0)
    return J
###############################################################################        
def polyactrow(u,w):
    
    dims_u = np.shape(u)
    u = np.reshape(u,(1,(dims_u[0]*dims_u[1])),order='f')
    u = np.transpose(u)
    h = np.zeros(dims_u[0],dims_u[1])
    h[u <= -w/2] = 0
    h[u >= w/2] = 1
    if w>0:
        h[(-w/2)<u] = 0.5*(((2/w)*u[-w/2<u]+1)**2)
        h[u<=0] = 0.5*(((2/w)*u[u<=0]+1)**2)
        h[0<u] = 1-0.5*(((2/w)*u[0<u]-1)**2)
        h[u<=(w/2)] = 1-0.5*(((2/w)*u[u<(w/2)]-1)**2)
    return h
###############################################################################        
def dpolyactrow(u,w):
    
    dims_u = np.shape(u)
    u = np.reshape(u,(1,(dims_u[0]*dims_u[1])),order='f')
    dh = np.zeros(dims_u[0],dims_u[1])
    dh[u <= -w/2] = 0
    dh[u >= w/2] = 1
    if w>0:
        dh[-w/2<u] = (4/(w**2))*u[-w/2<u]+(2/w)
        dh[u<=0] = (4/(w**2))*u[u<=0]+(2/w)
        dh[0<u] = -(4/(w**2))*u[0<u]+(2/w)
        dh[u<=(w/2)] = -(4/(w**2))*u[u<=(w/2)]+(2/w)
    return dh
            
        
        
        
        