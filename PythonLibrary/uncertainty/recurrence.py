import numpy as np
from scipy import special as sp

def jacobi_recurrence(N, alpha=0., beta=0., probability=False):
    # Returns the first N recurrence coefficient pairs for the Jacobi
    # polynomial family.

    if N < 1:
        return np.ones((0,2))

    ab = np.ones((N,2)) * np.array([beta**2.- alpha**2., 1.])

    # Special cases
    ab[0,0] = (beta - alpha) / (alpha + beta + 2.)
    ab[0,1] = np.exp( (alpha + beta + 1.) * np.log(2.) +
                      sp.gammaln(alpha + 1.) + sp.gammaln(beta + 1.) -
                      sp.gammaln(alpha + beta + 2.)
                    )

    if N > 1:
        ab[1,0] /= (2. + alpha + beta) * (4. + alpha + beta)
        ab[1,1] = 4. * (alpha + 1.) * (beta + 1.) / (
                   (alpha + beta + 2.)**2 * (alpha + beta + 3.) )

    inds = np.arange(2.,N)
    ab[2:,0] /= (2. * inds + alpha + beta) * (2 * inds + alpha + beta + 2.)
    ab[2:,1] = 4 * inds * (inds + alpha) * (inds + beta) * (inds + alpha + beta)
    ab[2:,1] /= (2. * inds + alpha + beta)**2 * (2. * inds + alpha + beta + 1.) * (2. * inds + alpha + beta - 1)

    ab[:,1] = np.sqrt(ab[:,1])

    if probability:
        ab[0,1] = 1.

    return ab

def hermite_recurrence(N, rho=0., probability=False):
    # Returns the first N recurrence coefficient pairs for the Hermite
    # polynomial family.

    if N < 1:
        return np.ones((0,2))

    ab = np.zeros((N,2))
    ab[0,1] = sp.gamma(rho+0.5)

    ab[1:,1] = 0.5*np.arange(1., N)
    ab[np.arange(N) % 2 == 1,1] += rho

    ab[:,1] = np.sqrt(ab[:,1])

    if probability:
        ab[0,1] = 1.

    return ab
