# Routines for Weighted Approximate Fekete Points

import numpy as np
from recurrence import jacobi_recurrence
from opolynd import opolynd_eval

def legendre_wafp(lambdas, M=1e3, sampler=None):
    """
    Generate M (= lambdas.shape[0]) weighted approximate Fekete points
    using randomized sampling.
    """

    from scipy.linalg import qr
    from legendre_induced import induced_distribution_mixture_sampling

    if lambdas.ndim == 1:
        lambdas = np.reshape(lambdas, [lambdas.size, 1])

    N, d = lambdas.shape

    ab = jacobi_recurrence(lambdas.max() + 1, alpha=0., beta=0., probability=True)

    if sampler is None:
        sampler = lambda MM: induced_distribution_mixture_sampling(lambdas, MM)

    # Choose at least 2*N samples
    M = max(M, 2*N)

    x = sampler(M)
    V = opolynd_eval(x, lambdas, ab)
    _, _, p = qr(V.T/np.sqrt(np.sum(V**2,axis=1)), pivoting=True, mode='economic')

    return x[p[:N], :]

def legendre_wafp_enrichment(x, lambdas, M_enrich, sampler=None):
    """
    Adds M_enrich points to the existing point set x by (approximate)
    determinant maximization.
    """

    from legendre_induced import induced_distribution_mixture_sampling

    if lambdas.ndim == 1:
        lambdas = np.reshape(lambdas, [lambdas.size, 1])

    if sampler is None:
        sampler = lambda MM: induced_distribution_mixture_sampling(lambdas, MM)

    M0 = x.shape[0]
    N, d = lambdas.shape
    ab = jacobi_recurrence(lambdas.max() + 1, alpha=0., beta=0., probability=True)

    M = max(M0, 2*N)

    while x.shape[0] < M0 + M_enrich:
        V = opolynd_eval(x, lambdas, ab)
        W = (V.T/np.sqrt(np.sum(V**2,axis=1))).T * np.sqrt(float(N)/float(x.shape[0]))
        G = np.dot(W.T, W)
        iG = np.linalg.inv(G)

        xs = sampler(M)
        Vs = opolynd_eval(xs, lambdas, ab)
        Ws = (Vs.T/np.sqrt(np.sum(Vs**2,axis=1))).T
        dets = np.sum((Ws*np.dot(Ws, iG))**2, axis=1)
        ind = np.argmax(dets)

        x = np.vstack([x, xs[ind,:]])

    return x

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    from indexing import total_degree_indices, hyperbolic_cross_indices
    from recurrence import jacobi_recurrence
    from opolynd import opolynd_eval

    d, k = 9, 7

    #lambdas = total_degree_indices(d, k)
    lambdas = hyperbolic_cross_indices(d, k)
    N = lambdas.shape[0]

    x = legendre_wafp(lambdas, M=3e3)

    ab = jacobi_recurrence(lambdas.max() + 1, alpha=0., beta=0., probability=True)

    V = opolynd_eval(x, lambdas, ab)
    W = (V.T/np.sqrt(np.sum(V**2,axis=1))).T * np.sqrt(float(N)/float(x.shape[0]))

    M_enrich = 20
    x2 = legendre_wafp_enrichment(x, lambdas, M_enrich, sampler=None)
    V2 = opolynd_eval(x2, lambdas, ab)
    W2 = (V2.T/np.sqrt(np.sum(V2**2,axis=1))).T * np.sqrt(float(N)/float(x2.shape[0]))

    # W: unenriched
    # W2: enriched
