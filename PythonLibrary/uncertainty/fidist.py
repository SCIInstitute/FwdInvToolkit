# Fast routines for induced distributions

import os, subprocess
import numpy as np
from scipy.io import loadmat

import opoly1d, recurrence

# Matlab binary location:
matlab_binary = '/Applications/MATLAB_R2015b.app/bin/matlab'
#matlab_binary = '/Applications/MATLAB_R2016a.app/bin/matlab'

# Abs path for data directory
cwd = os.path.dirname(os.path.abspath(__file__))
data_directory = os.path.join(cwd, '../data/')

def loadfile_jacobi(filename, alpha, beta, n):

    filename = os.path.join(data_directory, filename)
    try:
        data = loadmat(filename)['data'].flatten()
    except:
        data = np.zeros(0)

    if data.size < n+1:
        # Run matlab to generate/populate file

        print("Calling matlab....")
        cwd = os.path.dirname(os.path.abspath(__file__))

        command  = "cd(" + "'" + cwd + "'); cd ..; "

        command += "data = load_fjacobi({:d}, {:.4f}, {:.4f}); ".format(n, alpha, beta)
        command += "data = fidistinv_jacobi_setup({:d}, {:.4f}, {:.4f}, data); ".format(n, alpha, beta)
        command += "save_fjacobi(data, {:.4f}, {:.4f}); ".format(alpha, beta)

        command += "exit"
        #print command

        subprocess.call([matlab_binary, "-nodisplay", "-r", command])
        print("...finished")

        data = loadmat(filename)['data'].flatten()

    return data

def loadfile_hfreud(filename, alpha, rho, n):

    filename = os.path.join(data_directory, filename)
    try:
        data = loadmat(filename)['data'].flatten()
    except:
        data = np.zeros(0)

    if data.size < n+1:
        # Run matlab to generate/populate file

        print("Calling matlab....")
        cwd = os.path.dirname(os.path.abspath(__file__))

        command  = "cd(" + "'" + cwd + "'); cd ..; "

        command += "data = load_fhfreud({:d}, {:.4f}, {:.4f}); ".format(n, alpha, rho)
        command += "data = fidistinv_hfreud_setup({:d}, {:.4f}, {:.4f}, data); ".format(n, alpha, rho)
        command += "save_fhfreud(data, {:.4f}, {:.4f}); ".format(alpha, rho)

        command += "exit"
        print(command)

        subprocess.call([matlab_binary, "-nodisplay", "-r", command])
        print("...finished")

        data = loadmat(filename)['data'].flatten()

    return data


def fidistinv_jacobi(u, n, alpha, beta):

    assert alpha > -1. and beta > -1.
    filename = 'jacobi_{:3.4f}_{:3.4f}.mat'.format(alpha,beta)
    data = loadfile_jacobi(filename, alpha, beta, np.max(n))

    return fidistinv_driver(u, n, data)

def fidistinv_hfreud(u, n, alpha, rho):

    filename = 'half_freud_{:3.4f}_{:3.4f}.mat'.format(alpha,rho)
    data = loadfile_hfreud(filename, alpha, rho, np.max(n))

    return fidistinv_driver(u, n, data)

def fidistinv_freud(u, n, alpha, rho):

    x = np.zeros(u.shape)
    if n.size == 1:
        if n[0] % 2 == 0:
            x = np.sqrt( fidistinv_hfreud( np.abs(2.*u-1.), n/2,     alpha/2., (rho-1.)/2. ) )
        else:
            x = np.sqrt( fidistinv_hfreud( np.abs(2.*u-1.), (n-1)/2, alpha/2., (rho+1.)/2. ) )

    else:
        evenflags = n % 2 == 0
        oddflags = np.invert(evenflags)

        x[evenflags] = np.sqrt( fidistinv_hfreud( np.abs(2.*u[evenflags]-1.), n[evenflags]/2, \
                                                  alpha/2., (rho-1.)/2. ) )
        x[oddflags]  = np.sqrt( fidistinv_hfreud( np.abs(2.*u[oddflags]-1.), (n[oddflags]-1)/2, \
                                                  alpha/2., (rho+1.)/2. ) )

    # Reflect
    flags = u < 0.5
    x[flags] = -x[flags]
    return x

def fidistinv_driver(u, n, data):
    # Main driver for computing approximate (fast) induced
    # distributions.

    # Proximity to edge of interval at which we peg u to interval
    # endpoints
    inttol = 1e-8

    assert np.all(u>=0.) and np.all(u<=1.)
    assert np.all(n>=0)

    x = np.zeros(u.size)
    if u.size == 0:
        return x

    # Sort n (degree) values
    if n.size > 1:
        assert u.size == n.size
        inds = np.argsort(n)
        indsep = np.argwhere( np.diff(n[inds]) > 0 ) + 1
    else:
        inds = np.arange(0, u.size, dtype=int)
        indsep = np.zeros(0)

    # Create partition of indices in u
    indsep = np.insert(indsep, 0, 0)
    indsep = np.append(indsep, u.size)

    ushape = u.shape
    u = u.flatten()

    # Loop over each value of n (degree)
    for (i1,i2) in zip( indsep[:-1], indsep[1:] ):
        #nn = n[inds[i1]]
        #uu = u[inds[i1:i2]]

        x[inds[i1:i2]] = fidistinv_driver_helper(u[inds[i1:i2]], data[n[inds[i1]]])

    x.reshape(ushape)
    return x

def fidistinv_driver_helper(u, data):

    tol = 1e-12 # machine eps guarding


    # Construct Vandermonde-like matrix
    M = data.shape[0] - 6
    ab = recurrence.jacobi_recurrence(M+1, alpha=-1/2., beta=-1/2.)

    edges = np.concatenate([[-np.inf], data[0,:], [data[1,-1]], [np.inf]])
    bins = np.digitize(u, edges)
    B = edges.size-1;

    x = np.zeros(u.size)
    x[bins==1] = data[2,0]
    x[bins==B] = data[3,-1]

    for qb in range(2, B):
        mask = (bins==qb)
        if not(np.any(mask)):
            continue

        q = qb -2
        vgrid = (u[mask] - data[0,q])/(data[1,q] - data[0,q])*2. - 1.

        V = opoly1d.opoly1d_eval(vgrid, list(range(M)), ab)
        temp = np.dot(V, data[6:,q])
        temp /= ( (1. + vgrid)**(data[4,q]) * (1. - vgrid)**(data[5,q]) )

        if data[4,q] != 0.:
            # Set LHS to be 0
            flags = np.abs(u[mask] - data[0,q]) < tol
            temp[flags] = 0.
            temp *= (data[3,q] - data[2,q])
            temp += data[2,q]
        else:
            # Set RHS to be 0
            flags = np.abs(u[mask] - data[1,q]) < tol
            temp[flags] = 0.
            temp *= (data[3,q] - data[2,q])
            temp += data[3,q]

        x[mask] = temp

    return x

def idist_mixture_sampling_tensorial(M, d, k, disttype, **kwargs):
    # Returns M samples in d dimensions from a mixture of induced
    # distributions, where the mixture is associated to a degree-k
    # tensor-product distribution. The input disttype should be
    # "jacobi", "hfreud", or "freud". The keyword arguments required are
    # (alpha,beta) for jacobi, (alpha,rho) for hfreud, and (alpha,rho)
    # for freud.

    assert k >= 0   # non-negative degree
    assert d >= 1   # positive dimension

    # Tensorial is easy: all coordinates have the same distribution.
    u = np.random.rand(M*d)
    n = np.floor((k+1)*np.random.rand(M*d)).astype(int)
    n[n==(k+1)] = k

    if disttype == "jacobi":
        x = fidistinv_jacobi(u, n, **kwargs)
    elif disttype == "hfreud":
        x = fidistinv_hfreud(u, n, **kwargs)
    elif disttype == "freud":
        x = fidistinv_freud(u, n, **kwargs)
    else:
        assert False

    if d > 1:
        return x.reshape((M,d))
    else:
        return x

if __name__ == "__main__":

    import pdb
    from matplotlib import pyplot as plt
    import opolynd

    alpha = 0.
    beta = 0.
    N = 11
    u = np.linspace(0, 1, N)
    u[-1] -= 1e-6
    n = np.arange(N, dtype=int)

    # Testing loading jacobi data and evaluating fidistinv_jacobi
    filename = 'jacobi_{:3.4f}_{:3.4f}.mat'.format(alpha, beta)
    data = loadfile_jacobi(filename, alpha, beta, np.max(n))
    x = fidistinv_jacobi(u, n, alpha, beta)
    # Answer from matlab:
    xexact = [-1.000000000000000,  -0.928317766722557,  -0.887228979688486 , -0.517492607529973,  -0.210828858592612,   0.000000000000000, 0.377790521226704,   0.587803521443663,   0.842636353867568,   0.936846338676717,   1.000000000000000]
    print("Jacobi error: {:.4e}".format(np.linalg.norm(x - xexact)))

    # Testing creating other jacobi data with matlab
    alpha = np.random.rand(1)[0]*11 - 1.
    beta  = np.random.rand(1)[0]*11 - 1.
    filename = 'jacobi_{:3.4f}_{:3.4f}.mat'.format(alpha, beta)
    data = loadfile_jacobi(filename, alpha, beta, 2)

    # Testing half-line Freud evaluations
    alpha = 1.
    rho = 0.
    filename = 'half_freud_{:3.4f}_{:3.4f}.mat'.format(alpha, rho)
    data = loadfile_hfreud(filename, alpha, rho, np.max(n))
    x = fidistinv_hfreud(u, n, alpha, rho)
    # Answer from matlab:
    xexact = [ 0,   0.119552024763544,  1.230335315134000,   3.334709190393895,   6.400179446826883,  10.2875349395399, 18.458183625522008 , 23.995991765631384 , 29.564086690424123 , 35.880003534475335,  63.196084473861298]
    print("Laguerre error: {:.4e}".format(np.linalg.norm(x - xexact)))

    # Testing full-line Freud evaluations
    alpha = 2.
    rho = 0.
    x = fidistinv_freud(u, n, alpha, rho)
    # Answer from matlab:
    xexact = [ -5.042890705179727,  -1.523421749424712 , -1.668998069614480 , -1.728340302085805 , -0.993630459736906,                   0, 1.004265363591303,   2.206804803887973,   3.417302103995109,   4.003844362463767,   6.154799836463579]

    print("Hermite error: {:.4e}".format(np.linalg.norm(x - xexact)))

    # Testing random sampling from induced mixture (tensorial)
    M = 10000   # Number of samples
    d = 1       # Dimension
    k = 22      # Degree
    alpha = 0.  # Jacobi parameter
    beta = 0.   # Jacobi parameter
    x = idist_mixture_sampling_tensorial(M, d, k, "jacobi", alpha=alpha, beta=beta)

    plt.plot(np.sort(x), np.arange(x.size,dtype=float)/x.size)
    plt.title('CDF for Legendre, k={:d}'.format(k))

    k = 22
    alpha = 2.  # Freud (Hermite) parameter
    rho = 0.    # Freud (Hermite) parameter
    x = idist_mixture_sampling_tensorial(M, d, k, "freud", alpha=alpha, rho=rho)
    plt.figure()
    plt.plot(np.sort(x), np.arange(x.size,dtype=float)/x.size)
    plt.title('CDF for Hermite, k={:d}'.format(k))

    x = idist_mixture_sampling_tensorial(M, 5, k, "freud", alpha=alpha, rho=rho)

    # Testing that sampling like this is well-conditioned
    M = 2000
    k = 30
    d = 2
    disttype = "freud"
    alpha = 2.
    rho = 0.
    ab = recurrence.hermite_recurrence(k+1,rho=rho)

    # Sampler
    x = idist_mixture_sampling_tensorial(M, d, k, disttype, alpha=alpha, rho=rho)

    # Generate tensor-product multi-indices
    l = np.arange(k+1)
    lx,ly = np.meshgrid(l,l)
    L = np.concatenate((lx.reshape(lx.size,1), ly.reshape(ly.size,1)), axis=1)
    # Vandermonde-like matrix
    V = opolynd.opolynd_eval(x, L, ab)
    # Christoffel preconditioner
    V = (V.T* np.sqrt( (k+1.)**2. / np.sum(V**2, axis=1))).T

    plt.show()
