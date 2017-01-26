# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 11:11:58 2016

@author: andrewcorbato
"""

def linesearch(phi,dphidx,x,p,c1,c2,alphamin,alpharate,alphamax):
    # Algorithm: Line Search Algorithm (Algorithm 3.2, pg.59 of "Numerical
    # Optimization" by Nocedal & Wright
    # Author: Burak Erem
    # ----------------------------------------
    # INPUT: phi = objective function
    #        dphidx = gradient of objective function
    #        p = search direction
    #        c1 = "sufficient decrease condition" parameter
    #        c2 = "curvature condition" parameter
    #        alphamin = minimum starting alpha (should be >0)
    #        alpharate = multiplicative rate at which alpha grows
    #                    (i.e. next=alpharate*prev)
    #        alphamax = upper bound on the chosen alpha
    # ----------------------------------------
    # OUTPUT: alphafinal = chosen alpha
    import numpy as np
    import scipy as sp
    phialpha = lambda a: phi[x+a*p]
    dphidalpha = lambda a: np.transpose(dphidx[x+a*p])*p
    
    # "Exact" line search
    # alphas=logspace(log10(alphamin),log10(alphamax),5);
    # alphas=linspace(alphamin,alphamax,ceil(numel(x)/10));
    # need to update comments
    
    eps = np.finfo(float).eps
    alphas = np.linspace(eps,alphamax,10)
    dims_alphas = np.shape(alphas)
    costs = np.zeros(dims_alphas[0]*dims_alphas[1])
    numel_alphas = dims_alphas[0]*dims_alphas[1]
    
    for ii in range(0,numel_alphas):
        costs[ii] = phialpha(alphas[ii])
    ind = np.argmin(costs,axis=0)
    alphafinal = alphas[ind]
    return alphafinal
    
    # initialize the algorithm
    alpha[0] = 0
    alpha[1] = alphamin
    i = 1
    while True:
        # check for the first terminating condition
        currphi = phialpha(alpha[i])
        if (currphi>(phialpha(0)+c1*alpha[i]*dphialpha[0]) or currphi>=phialpha(alpha(i-1)) and i>1):
            alphafinal = zoom(alpha[i-1],alpha[i],phialpha,dphidalpha,c1,c2)
            # return
        
        # check for the second condition
        currdphialpha = dphidalpha(alpha[i])
        if np.abs(currdphialpha) <= -c2*dphidalpha(0):
            alphafinal = alpha[i]
            # return
        
        # check for the third terminating condition
        if currdphidalpha>=0:
            alphafinal = zoom(alpha[i],alpha[i-1],phialpha,dphidalpha,c1,c2)
            # return;
        
        # if no termination yet, update alpha
        alpha[i+1] = alpharate*alpha[i]
        
        # check for last terminating condition
        if alpha[i+1] > alphamax:
            alphafinal = alphamax
            # return
        i += 1
###############################################################################    
# this version of zoom uses cubic interpolation
def zoom(alphalo,alphahi,phi,dphi,c1,c2):
    
    import numpy as np
    
    while True:
        # cubic interpolation
        d1 = dphi[alphalo]+dphi[alphahi]-3*(phi[alphalo]-phi[alphahi])/(alphalo-alphahi)
        d2 = np.sqrt(np.max(d1**2-dphi[alphalo]*dphi[alphahi,0]))
        alpha = alphalo-(alphalo-alphahi)*(dphi[alphahi]+d2-d1)/(dphi[alphahi]-dphi[alphalo]+2*d2)
        
        # sanity checks on value of alph
        if alpha>alphahi:
            alpha = alphahi
            # return
        if alpha<alphalo:
            alpha = alphalo
            # return
        if np.abs(alphahi-aplhalo)<= eps:
            alpha = alphalo
            # return
            
        # remainder of zoom algorithm
        currphi = phi[alpha]
        if (currphi>(phi[0]+c1*alpha*dphi[0]) or currphi >= phi(alphalo)):
            alphahi = alpha
        else:
            currdphi = dphi[alpha]
            if np.abs*currdphi <= -c28dphi[0]:
                # return
            if currdphi*(alphahi-alphalo)>= 0:
                alphahi = alphalo
                
            alphalo = alpha
            
return alpha
