import numpy as np

def gauss_quadrature(ab, N):
    """
    Computes the N-point Gauss quadrature rule associated to the
    recurrence coefficients ab.
    """

    from numpy.linalg import eigh

    n = ab.shape[0]
    if n+1 < N:
        raise IndexError('Require N+1 recurrence coefficients for an N-point rule.')

    J = np.diag(ab[1:N,1], k=1) + np.diag(ab[:N,0],k=0) + np.diag(ab[1:N,1], k=-1)
    lamb,v = eigh(J)

    return lamb, v[0,:]**2

if __name__ == "__main__":

    from recurrence import jacobi_recurrence
    from opoly1d import opoly1d_eval

    alpha = 0
    beta = 0

    N = 12
    ab = jacobi_recurrence(N+1, alpha=alpha, beta=beta, probability=True)

    x,w = gauss_quadrature(ab,N)

    V = opoly1d_eval(x, list(range(N)), ab)

    G = np.dot(np.dot(V.T, np.diag(w)), V)
    err = np.linalg.norm(G - np.eye(G.shape[0]))
