import os.path
import numpy as np

def load_sparse_grid(dim, level):
    """
    Loads sparse grid data from saved ascii files. Relative to this
    file's path, searches in a subdirectory "sgdata" for a file whose
    name is

      'dim' + str(dim) + '_l' + str(level) + '.dat'

    This file should be plaintext file with dim+1 columns. The first
    column contains quadrature weights and the remaining dim columns
    contain nodal coordinates.
    """

    relative_path = os.path.dirname(os.path.realpath(__file__))
    filename = 'dim' + str(dim) + '_l' + str(level) + '.dat'
    filename = os.path.join(relative_path,  'sgdata', filename)

    try:
        wx = np.loadtxt(filename)
    except IOError:
        print(('Cannot load data: file ' + filename + ' does not exist'))
        return

    return wx[:,1:], wx[:,0]

if __name__ == "__main__":

    import sensitivity, indexing, test_functions, opoly1d, opolynd

    dim = 3
    level = 3;

    x,w = load_sparse_grid(dim, level)

    # Test on a multidimensional genz function
    c = np.random.randn(dim)  # Frequency parameter
    r = np.random.randn(1)    # phase pharameter

    f = test_functions.genz_oscillatory(c, r)
    fx = f(x)

    # Compute integral of scalar function
    Q_sparse_grid = np.sum(w*fx)
    Q_exact = test_functions.genz_oscillatory_integral(c, r)
    Q_relerror = np.abs(Q_sparse_grid - Q_exact)/np.abs(Q_exact)

    print(("Sparse grid integration error for Genz oscillatory test function is {:1.4e} using a level {:d} grid in {:d} dimensions".format(Q_relerror[0], level, dim)))

    # Compute a polynomial chaos expansion (PCE) assuming a uniform
    # distribution

    order = level

    # Unifrom distribution
    ab = opoly1d.jacobi_recurrence(order+1, alpha=0., beta = 0., probability=True)
    lambdas = indexing.total_degree_indices(dim, order)
    V = opolynd.opolynd_eval(x, lambdas, ab)

    # Test integration accuracy of sparse grid:
    G = np.dot(np.dot(np.transpose(V), np.diag(w)), V)
    print(("Sparse grid integration error for orthogonal polynomials up to degree {:d} is {:1.4e} using a level {:d} grid in {:d} dimensions".format(order, np.linalg.norm(G - np.eye(G.shape[0])), level, dim)))

    # Compute PCE coefficients
    coeffs = np.dot(np.transpose(V), w*fx)

    # Now things like sensitivity can be computed:

    # Total sensitivity for each dimension:
    pce_TS = sensitivity.pce_total_sensitivity(coeffs, lambdas, list(range(dim)))
    # If you compare pce_TS with c, the qualitative trend in these
    # vectors is similar

    # Global sensitivity for dimension 0, for dimensions [0,1], and for
    # dimensions [0,2]
    pce_GS = sensitivity.pce_global_sensitivity(coeffs, lambdas, [[0,], [0,1], [0,2]])

    # All the global sensitivities in 3 dimensions
    if dim == 3:
        dim3_power_set = [[0,], [1,], [2,], [0,1], [1,2], [0,2], [0,1,2]]
        pce_all_GS = sensitivity.pce_global_sensitivity(coeffs, lambdas, dim3_power_set)
        # The global sensitivities satisfy np.sum(pce_all_GS) == 1.
