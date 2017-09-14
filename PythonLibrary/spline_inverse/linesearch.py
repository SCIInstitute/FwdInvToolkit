import numpy as np
import scipy as sp


def linesearch(phi,dphidx,x,p,alphamax):
    # Algorithm: Line Search Algorithm (Algorithm 3.2, pg.59 of "Numerical
    # Optimization" by Nocedal & Wright
    # Author: Burak Erem
    # ----------------------------------------
    # INPUT: phi = objective function
    #        dphidx = gradient of objective function
    #        p = search direction
    #        alphamax = upper bound on the chosen alpha
    # ----------------------------------------
    # OUTPUT: alphafinal = chosen alpha

    phialpha = lambda a: phi(x+a*p)
    dphidalpha = lambda a: dphidx(x+a*p).T @ p
    
    # "Exact" line search    
    eps = np.finfo(float).eps
    alphas = np.linspace(eps,alphamax,10)

    costs = np.zeros(alphas.shape)    
    for i in range(0,alphas.size):
        costs[i] = phialpha(alphas[i])

    ind = np.argmin(costs, axis=0)
    alphafinal = alphas[ind]

    return alphafinal
    
if __name__ == '__main__':
    N = 2
    phi = lambda x: np.sum(x**2)
    dphidx = lambda x: 2*np.ones((N,1))
    x = np.ones((N,1))
    p = -np.ones((N,1))
    alphamax = 10

    alphafinal = linesearch(phi,dphidx,x,p,10)
    print(alphafinal)

