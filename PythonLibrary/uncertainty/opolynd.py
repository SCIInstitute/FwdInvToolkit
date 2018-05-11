import numpy as np
import opoly1d

def opolynd_eval(x, lambdas, ab):
    # Evaluates tensorial orthonormal polynomials associated with the
    # univariate recurrence coefficients ab.

    try:
        M, d = x.shape
    except Exception:
        d = x.size
        M = 1
        x = np.reshape(x, (M, d))

    N, d2 = lambdas.shape

    assert d==d2, "Dimension 1 of x and lambdas must be equal"

    p = np.ones([M, N])

    for qd in range(d):
        p = p * opoly1d.opoly1d_eval(x[:,qd], lambdas[:,qd], ab)

    return p

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D
    import opoly1d, indexing

    d = 2
    k = 4

    ab = opoly1d.jacobi_recurrence(k+1, probability=True)

    N = 50
    x = np.linspace(-1, 1, N)
    X,Y = np.meshgrid(x,x)

    XX = np.concatenate((X.reshape(X.size,1), Y.reshape(Y.size,1)), axis=1)

    lambdas = indexing.total_degree_indices(d, k)

    p = opolynd_eval(XX, lambdas, ab)

    j = 11
    assert j < lambdas.shape[0]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, p[:,j].reshape(N,N), cmap=cm.coolwarm, linewidth=0,antialiased=True)
    plt.show()
