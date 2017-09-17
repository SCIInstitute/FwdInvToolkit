import numpy as np
from linesearch import linesearch


def steepestdescent(phi, dphidx, xinit, mingradnorm, alphamax=10):
    # Algorithm: Steepest Descent
    # Author: Burak Erem
    # --------------------------------------
    # INPUT: phi = objective function
    #        dphidx = gradient of objective function
    #        xinit = initialization of algorithm step sequence
    #        mingradnorm = smallest 2-norm of dphidx before stopping
    # --------------------------------------
    # OUTPUT: x = matrix whose columns are the sequence of x points
    #         dfdx = matrix whose columns are the gradients at x
    #         alpha = array of the step size parameters used
    #

    # Initialization:
    k = 0
    x = []
    x.append(xinit)
    f = phi(x[k])
    dfdx = []
    dfdx.append(dphidx(x[k]))
    p = []
    p.append(-dfdx[k])

    alpha = []

    costs = []
    costs.append(phi(x[k]))
    cost_reduction = np.inf

    print('Init obj fun={0}'.format(costs[k]))

    # The main loop
    while (cost_reduction > mingradnorm):
        # Line Search (satisfies strong Wolfe conditions)
        alpha.append(linesearch(phi, dphidx, x[k], p[k], alphamax))

        # Set next sequence step x_(k+1)
        x.append(x[k] + alpha[k] * p[k])

        # Evaluate gradient(f_(k+1))
        dfdx.append(dphidx(x[k + 1]))

        # Set search direction p_(k+1)
        p.append(-dfdx[k + 1])
        costs.append(phi(x[k + 1]))
        cost_reduction = costs[k] - costs[k + 1]

        # Report progress on the current iteration
        print('k={0}    objfun={1}      costreduction={2}'.format(k, costs[k + 1], cost_reduction))

        # Increment k
        k = k + 1

    return x, dfdx, alpha


if __name__ == '__main__':
    N = 2
    phi = lambda x: np.sum(x ** 2)
    dphidx = lambda x: 2 * np.ones((N, 1))
    xinit = np.ones((N, 1))
    alphamax = 0.1
    mingradnorm = 1e-15

    x, dfdx, alpha = steepestdescent(phi, dphidx, xinit, mingradnorm, alphamax)
    print(x[-1])


