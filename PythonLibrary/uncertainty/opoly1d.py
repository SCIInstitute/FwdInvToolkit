import numpy as np
from scipy import special as sp

def opoly1d_eval(x, n, ab, jac=0):
    # Evaluates univariate orthonormal polynomials given their
    # three-term recurrence coefficients ab.
    #
    # Evaluates the jac'th derivative. (Default = 0)
    #
    # Returns a numel(x) x numel(n) x numel(jac) array.

    n = np.asarray(n)
    if n.size < 1 or x.size < 1:
        return np.zeros(0)

    nmax = np.max(n)
    assert nmax < ab.shape[0]
    assert np.min(n) > -1
    assert jac >= 0

    p = np.zeros( x.shape + (nmax+1,) )
    xf = x.flatten()

    p[:,0] = 1/ab[0,1]

    if nmax > 0:
        p[:,1] = 1/ab[1,1] * ( (xf - ab[0,0])*p[:,0] )

    for j in range(2, nmax+1):
        p[:,j] = 1/ab[j,1] * ( (xf - ab[j-1,0])*p[:,j-1] - ab[j-1,1]*p[:,j-2] )

    p = p[:,n.flatten()]

    if type(jac) == int:
        if jac == 0:
            return p
        else:
            jac = [jac]

    preturn = np.zeros([p.shape[0], p.shape[1], len(jac)])
    for (qi,qval) in enumerate(jac):
        if qval == 0:
            preturn[:,:,qi] = p

    for qd in range(1, max(jac)+1):
        pd = np.zeros(p.shape)

        for qn in range(qd,nmax+1):
            if qn == qd:
                # The following is an over/underflow-resistant way to
                # compute ( qd! * kappa_{qd} ), where qd is the
                # derivative order and kappa_{qd} is the leading-order
                # coefficient of the degree-qd orthogonal polynomial.
                # The explicit formula for the lading coefficient of the
                # degree-qd orthonormal polynomial is prod(1/b[j]) for
                # j=0...qd.
                pd[:,qn] = np.exp( sp.gammaln(qd+1) - np.sum( np.log( ab[:(qd+1),1] ) ) )
            else:
                pd[:,qn] = 1/ab[qn,1] * ( ( xf - ab[qn-1,0] ) * pd[:,qn-1] - ab[qn-1,1] * pd[:,qn-2] )
        for (qi,qval) in enumerate(jac):
            if qval == 0:
                preturn[:,:,qi] = pd

        p = pd

    return preturn.flatten()

def christoffel_function(x, k, ab):
    # Computes the normalized (inverse) Christoffel function lambda, defined as
    #
    #   lambda**2 = k / sum(p**2, axi=1),
    #
    # where p is a matrix containing evaluations of an orthonormal
    # polynomial family up to degree k-1, defined by the recurrence
    # coefficients ab.

    assert k > 0

    p = opoly1d_eval(x, list(range(k)), ab)
    return np.sqrt(float(k) / np.sum(p**2, axis=1))

def qpoly1d_eval(x, n, ab, jac=0):
    # Evalutes Christoffel-function normalized polynomials. These are
    # given by
    #
    #   q_k(x) = p_k(x) / sqrt( sum_{j=0}^{n-1} p_j^2 ), k = 0, ..., n-1
    #
    # The output is a x.size x n array

    assert n > 0

    q = np.zeros((x.size, n))
    q[:,0] = 1.
    qt = np.zeros(x.size)

    if n > 1:
        qt = 1/ab[1,1] * (x - ab[0,0]) * q[:,0]
        q[:,1] = qt / np.sqrt(1 + qt**2)

    for j in range(1, n-1):
        qt = 1/ab[j+1,1] * ( (x - ab[j,0])*q[:,j] - ab[j,1] * q[:,j-1] / np.sqrt(1 + qt**2) )
        q[:,j+1] = qt / np.sqrt(1 + qt**2)

    if type(jac) == int:
        if jac == 0:
            return q
        else:
            jac = [jac]

    assert False
    qreturn = np.zeros([q.shape[0], q.shape[1], len(jac)])
    for (qi,qval) in enumerate(jac):
        if qval == 0:
            qreturn[:,:,qi] = q

    for qd in range(1, max(jac)+1):
        asdf

    return qreturn

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    from recurrence import jacobi_recurrence

    alpha = 0.
    beta = np.pi
    nmax = 35
    probability_measure = True

    ab = jacobi_recurrence(nmax+1,alpha=alpha,beta=beta,probability=probability_measure)
    x = np.linspace(-1, 1, 300)

    p = opoly1d_eval(x, np.arange(nmax), ab)
    plt.plot(x, p[:,:10])
    plt.title('Polynomials')

    plt.figure()
    q = qpoly1d_eval(x, nmax, ab)
    plt.plot(x, q[:,:10])
    plt.title('Q polynomials')

    err = np.linalg.norm( q - p / np.sqrt( np.cumsum( p**2, axis=1 ) ) )

    plt.show()
