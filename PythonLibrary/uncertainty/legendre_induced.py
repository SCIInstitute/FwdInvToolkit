# Computes quantities assocated with sampling from an induced
# distribution mixture for the Legendre polynomials

import numpy as np

def induced_distribution_quadrule(n):
    """
    Returns a Gaussian quadrature rule that would be used by
    induced_distribution.
    """

    from recurrence import jacobi_recurrence
    from quadrature import gauss_quadrature

    n = np.atleast_1d(n)
    N = n.max()
    ab = jacobi_recurrence(N+2, alpha=0., beta=0., probability=True)
    return gauss_quadrature(ab, N+1)

def induced_distribution(x, n, quadrule=None):
    """
    Computes the induced distribution function, which is

    F_n(x) = 1/2 * \int_{-1}^x p_n^2(s) ds,

    where p_n is the degree-n orthonormal Legendre polynomial. With this
    normalization, the endpoint values are F_n(-1) = 0, and F_n(1) = 1
    for all n.

    This function is vectorized in both x and n, and supports the
    calling syntaxes

        induced_distribution(array, array)
        induced_distribution(array, scalar)

    In the first syntax, both arrays must have the same number of
    elements.

    If specified, the optional input quadrule is a list, [x,w],
    containing a quadrature rule that is assumed accurate for all
    polynomials up to degree 2*np.max(n.flatten()). If not given, this
    quadrature rule is generated as an appropriate-sized Gaussian
    quadrature rule.
    """

    from recurrence import jacobi_recurrence
    from opoly1d import opoly1d_eval

    x = np.atleast_1d(x)
    n = np.atleast_1d(n)
    if np.size(x) != np.size(n):
        if n.size != 1:
            raise ValueError("If n is the not the same size as x, then n must be a scalar.")

    assert (n >= 0).all(), "Degrees n must be non-negative"

    N = n.max()

    if quadrule is not None:
        xg, wg = quadrule[0], quadrule[1]
    else:
        xg, wg = induced_distribution_quadrule(n)

    # Output array
    u = np.zeros(x.shape)

    ab = jacobi_recurrence(N+2, alpha=0., beta=0., probability=True)

    # Unfortunately, I don't see an efficient way to avoid a for loop here
    if n.size > 1:
        for ind,xv in np.ndenumerate(np.asarray(x)):
            if xv <= -1:
                u[ind] = 0.
            elif xv >= 1:
                u[ind] = 1.
            else:
                vg = (xg+1.)/2.*(xv+1) - 1.
                u[ind] = (xv+1.)/2*np.dot((opoly1d_eval(vg, n[ind], ab)**2).T, wg)
    else:
        for ind,xv in np.ndenumerate(np.asarray(x)):
            if xv <= -1:
                u[ind] = 0.
            elif xv >= 1:
                u[ind] = 1.
            else:
                vg = (xg+1.)/2.*(xv+1) - 1.
                u[ind] = (xv+1.)/2*np.dot((opoly1d_eval(vg, n, ab)**2).T, wg)

    return u

def induced_distribution_inverse(u, n):
    """
    Computes the inverse of the induced distribution function. Any
    elements of u that are less than 0 get mapped to -1, and elements
    bigger than 1 get mapped to 1. (This is essentially a silent error
    correction.)

    Supports the same vectorization syntax as induced_distribution.
    """

    from scipy.optimize import bisect
    optoptions = {'xtol':1e-10, 'rtol':1e-10}

    u = np.atleast_1d(u)
    n = np.atleast_1d(n)

    x = np.zeros(u.shape)

    qrule = induced_distribution_quadrule(n)

    if n.size > 1:
        for ind,uv in np.ndenumerate(u):
            if uv <= 0:
                x[ind] = -1.
            elif uv >= 1:
                x[ind] = 1.
            else:
                f = lambda xx: induced_distribution(xx, n[ind], qrule) - uv
                x[ind] = bisect(f, -1., 1., **optoptions)
    else:
        for ind,uv in np.ndenumerate(u):
            if uv <= 0:
                x[ind] = -1.
            elif uv >= 1:
                x[ind] = 1.
            else:
                f = lambda xx: induced_distribution(xx, n, qrule) - uv
                x[ind] = bisect(f, -1., 1., **optoptions)

    return x

def induced_distribution_mixture_sampling(lambdas, M):
    """
    Returns M d-variate samples from the order-lambdas induced
    distribution. Each row of lambdas is a d-dimensional multi-index.
    """

    if lambdas.ndim == 1:
        lambdas = np.reshape(lambdas, [lambdas.size, 1])

    d = lambdas.shape[1]
    u = np.random.uniform(0., 1., [M,d])

    # Randomly sample M multi-indices
    lambdas_sampled = lambdas[np.random.randint(0,lambdas.shape[0],M),:]

    return induced_distribution_inverse(u, lambdas_sampled)


if __name__ == "__main__":

    from matplotlib import pyplot as plt

    import pdb
    from recurrence import jacobi_recurrence
    from opoly1d import opoly1d_eval
    import indexing

    n = 3
    x = np.linspace(-1., 1., 1e2)

    u = induced_distribution(x, n)

    qrule = induced_distribution_quadrule(n)
    u2 = induced_distribution(x, n, quadrule=qrule)
    uerr = np.linalg.norm(u-u2)

    x2 = induced_distribution_inverse(u, n)
    xerr = np.linalg.norm(x-x2)

    # Now do multiple n's
    N = 8
    ns = np.zeros(1000, dtype=int)
    ns = np.around(8*np.random.uniform(0, 1, 1000)).astype(int)
    xs = np.random.uniform(-1, 1, 1000)
    us = induced_distribution(xs, ns)
    xs2 = induced_distribution_inverse(us, ns)
    xserr = np.linalg.norm(xs-xs2)

    plt.plot(x, u, 'r', linewidth=3)
    plt.xlabel(r'$x$')

    plt.ylabel(r'$u = F_n(x)$')
    plt.show()

    d, k = 2, 5
    lambdas = indexing.total_degree_indices(d, k)
    x = induced_distribution_mixture_sampling(lambdas, 1e3)
