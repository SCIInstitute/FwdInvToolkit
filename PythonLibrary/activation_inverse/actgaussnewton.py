# -*- coding: utf-8 -*-
import scipy as sp
import numpy as np
import scipy.io as sio

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

#   import scipy and numpy package
#    import scipy as sp
#    import numpy as np

    #from sp.sparse import csr_matrix
    #L = csc_matrix.todense(L)
    
    
    tau = np.array(tauinit)
    A = np.array(A)
    Y = np.array(Y)
    L = np.array(L)
    
    dims_tau = np.shape(tau)
    step = 2*np.ones(dims_tau)
    print('size step = ',step.shape)
#   num = linesteps the number of steps in [0,1] to consider for the line search
    alpha = np.linspace(0,1, num = 100)
    dims_alpha = np.shape(alpha)
    dims_N = np.shape(A)
    N = dims_N[1]
    dims_T = np.shape(Y)
    T = dims_T[1]
    u = np.linspace(1,T,T)

#   if width variable isn't an array, make it a constant array
# need reshape
    #dims_w = len(w)
    #print(dims_w)
    #if dims_w == 1:
    w = w * np.ones((N,1))
    Iter = 0
        
    while np.linalg.norm(step)> minstep:
        print('starting while loop')
        
#       calculate the jacobian matrix and the residuals
        J = agnjacobian(A,Y,L,tau,Lambda,w)
        r = agnresidual(A,Y,L,tau,Lambda,w)
    
#       calculate the step direction
        G = np.dot(np.transpose(J),J)
        #print('trying to save')
        #sio.savemat('/Users/jess/Downloads/G.mat', {'G_tmp': G})
        #sio.savemat('/Users/jess/Downloads/r.mat', {'r_tmp': r})
        print('size r = ',r.shape)
        print('norm r = ',np.linalg.norm(r))
        print('max r = ',np.max(r))
        
        size_G = np.shape(G)
        print('size G = ',size_G)
        #print('G= ',G)
        print('norm G = ',np.linalg.norm(G))
        print('max G = ',np.max(G))
        s = np.linalg.svd(G, full_matrices = 0,compute_uv=0)
        cond_g = s[0]/s[len(s)-1]
        print('rcond of G =',1/cond_g)
        # improve the conditioning of matrix G
        cnt = 0
        print('starting condition test')
        while (1/np.linalg.cond(G,1))<=(np.finfo(float).eps):
            cnt = cnt +1
            G = G+Lambda*np.identity(size_G)
            print(cnt)

        print('conditioning of G complete')
        
        

#step = np.dot(np.linalg.lstsq(-G,np.transpose(J))[0],np.transpose(r))
        #find least sqares solution with qr decomp.
        qdc,rdc = np.linalg.qr(-G)
        step = np.dot(np.dot(np.linalg.inv(rdc),np.dot(qdc.T,J.T)),r.T)
        step = np.reshape(step,dims_tau)
        print('size step = ',step.shape)
        #sio.savemat('/Users/jess/Downloads/step.mat', {'step_tmp': step})

        print('step computed')

#       perform a line search in the step direction
        err = np.zeros(dims_alpha)
        for k in range(0,dims_alpha[0]):
            H = np.zeros((N,T))
            for n in range(0,N):
                H[n,:] = polyactrow(u-alpha[k]*step[n]-tau[n],w[n])
                
            err[k] = (np.linalg.norm(Y-np.dot(A,H),'fro'))**2+Lambda*(np.linalg.norm(np.dot(L,H),'fro'))**2
            if k>0:
                if err[k]>err[k-1]:
                    err = err[0:k]
                    break;

        print('for loop comlete')


           
#       update step with the result of the line search
        dims_err = np.shape(err)
        vec_err = err.T
        ind = np.argmin(err)
        print('ind = ',ind)
        step = alpha[ind]*step
        print('size step = ',step.shape)
        print('size tau = ',tau.shape)
        #        print('tau = ',tau)
#       update tau with the result of the step
        tau = tau+step
        print('size tau = ',tau.shape)
        #        print('tau = ',tau)
        
#       display the progress of the optimization routine
        Iter = Iter+1
        print('Iter = ',Iter)
        print('Step: %i\nSize: %f'%(Iter,np.linalg.norm(step)))
        return tau


###############################################################################
def agnresidual(A,Y,L,tau,Lambda,w):

    dim_A = np.shape(A)
    dim_Y = np.shape(Y)
    N = dim_A[1]
    T = dim_Y[1]
    
    u = np.linspace(1,T,T)
    
    H = np.zeros((N,T))
    
    for k in range(0,N):
      H[k,:] = polyactrow(u-tau[k],w[k])
    
    #sio.savemat('/Users/jess/Downloads/H.mat', {'H_tmp': H})

#print("size Y = ",np.shape(Y))
#    print("size A = ",np.shape(A))
#    print("size H = ",np.shape(H))

    E = Y-np.dot(A,H)
    R = np.sqrt(Lambda)*np.dot(L,H)
    dims_E = np.shape(E)
    dims_R = np.shape(R)
    vec_E = np.transpose(np.reshape(np.transpose(E),(dims_E[0]*dims_E[1])))
    vec_R = np.transpose(np.reshape(np.transpose(R),(dims_R[0]*dims_R[1])))

#    sio.savemat('/Users/jess/Downloads/E.mat', {'E_tmp': E})
#    sio.savemat('/Users/jess/Downloads/R.mat', {'R_tmp': R})
#    sio.savemat('/Users/jess/Downloads/vec_E.mat', {'vec_E_tmp': vec_E})
#    sio.savemat('/Users/jess/Downloads/vec_R.mat', {'vec_R_tmp': vec_R})

    r = np.concatenate((vec_E,vec_R), axis=0)
    return r

                 
###############################################################################
def agnjacobian(A,Y,L,tau,Lambda,w):
    dim_A = np.shape(A)
    dim_Y = np.shape(Y)
    M = dim_A[0]
    N = dim_A[1]
    T = dim_Y[1]
    dE = np.empty((M,T,N))
    dR = np.empty((N,T,N))
    HdH = np.zeros((1,T,N))
    u = np.linspace(1,T,T)
    
    for k in range(0,N):
      HdH[0,:,k] = dpolyactrow(u-tau[k],w[k])
    
    #    sio.savemat('/Users/jess/Downloads/HdH.mat', {'HdH_tmp': HdH})
    
    #print('size HdH = ',np.shape(HdH))
    #print('norm HdH = ',np.linalg.norm(HdH))
    #print('max HdH = ',np.max(HdH))

    A = np.reshape(A,(M,1,N)) # MxTxL
    L = np.reshape(L,(N,1,N)) # NxTxL
    
    HdHm = np.tile(HdH,[M,1,1]) # MxTxL
    HdHn = np.tile(HdH,[N,1,1]) # NxTxL

#    sio.savemat('/Users/jess/Downloads/HdHm.mat', {'HdHm_tmp': HdHm})
#    sio.savemat('/Users/jess/Downloads/HdHn.mat', {'HdHn_tmp': HdHn})

    A = np.tile(A, [1,T,1])
    L = np.tile(L, [1,T,1])

#sio.savemat('/Users/jess/Downloads/A.mat', {'A_tmp': A})
#   sio.savemat('/Users/jess/Downloads/L.mat', {'L_tmp': L})

    print("size A = ",np.shape(A))
#print('norm A = ',np.linalg.norm(A))
#    print('max A = ',np.max(A))
#    print("size L = ",np.shape(L))
#    print('norm L = ',np.linalg.norm(L))
#    print('max L = ',np.max(L))
#    print("size HdHm = ",np.shape(HdHm))
#    print('norm HdHm = ',np.linalg.norm(HdHm))
#    print('max HdHm = ',np.max(HdHm))
#    print("size HdHn = ",np.shape(HdHn))
#    print('norm HdHn = ',np.linalg.norm(HdHn))
#    print('max HdHn = ',np.max(HdHn))

    dE = np.multiply(A,HdHm)
    dR = -np.sqrt(Lambda)*np.multiply(L,HdHn)

#sio.savemat('/Users/jess/Downloads/dE.mat', {'dE_tmp': dE})
#    sio.savemat('/Users/jess/Downloads/dR.mat', {'dR_tmp': dR})
    
    print("size dE = ",np.shape(dE))
#    print('norm dE = ',np.linalg.norm(dE))
#    print('max dE = ',np.max(dE))
#    print("size dR = ",np.shape(dR))
#    print('norm dR = ',np.linalg.norm(dR))
#    print('max dR = ',np.max(dR))

    J = np.zeros(((M+N)*T,N))
    for k in range(0,N):
        tempE = dE[:,:,k]
        tempR = dR[:,:,k]
        dims_tempE = np.shape(tempE)
        dims_tempR = np.shape(tempR)
        vec_tempE = np.transpose(np.reshape(np.transpose(tempE),(dims_tempE[0]*dims_tempE[1])))
        vec_tempR = np.transpose(np.reshape(np.transpose(tempR),(dims_tempR[0]*dims_tempR[1])))
        #sio.savemat('/Users/jess/Downloads/tempE.mat', {'tempE_tmp': tempE})
        #sio.savemat('/Users/jess/Downloads/tempR.mat', {'tempR_tmp': tempR})
        #sio.savemat('/Users/jess/Downloads/vec_tempE.mat', {'vec_tempE_tmp': vec_tempE})
        #sio.savemat('/Users/jess/Downloads/vec_tempR.mat', {'vec_tempR_tmp': vec_tempR})
        J[:,k] = np.concatenate((vec_tempE,vec_tempR), axis=0)
        #sio.savemat('/Users/jess/Downloads/J1.mat', {'J_tmp': J})
        #return false
    #print('norm J = ',np.linalg.norm(J))
    #print('max J = ',np.max(J))
    #print('size J = ',np.shape(J))
    #print('J= ',J)
#sio.savemat('/Users/jess/Downloads/J.mat', {'J_tmp': J})
    return J
###############################################################################        
def polyactrow(u,w):
    
    dims_u = np.shape(u)
    h = np.zeros(dims_u)
    
    h[u <= -w/2] = 0
    h[u >= w/2] = 1
    if w>0:
        h[((-w/2)<u) & (u<=0)] = 0.5*(((2/w)*u[((-w/2)<u) & (u<=0)]+1)**2)
        h[(0<u) & (u<=(w/2))] = 1-0.5*(((2/w)*u[(0<u) & (u<=(w/2))]-1)**2)
    return h
###############################################################################        
def dpolyactrow(u,w):
  
    dims_u = np.shape(u)
    dh = np.zeros(dims_u)

#print('size u = ',dims_u)
#    print('u = ',u)


    dh[u <= -w/2] = 0
    dh[u >= w/2] = 1
    #    print('dh = ', dh)

    if w>0:
        dh[(-w/2<u) & (u<=0)] = (4/(w**2))*u[(-w/2<u) & (u<=0)]+(2/w)
        dh[(0<u) & (u<=(w/2))] = -(4/(w**2))*u[(0<u) & (u<=(w/2))]+(2/w)
          #print('dh = ', dh)
#return false
    return dh
            
        
        
        
        