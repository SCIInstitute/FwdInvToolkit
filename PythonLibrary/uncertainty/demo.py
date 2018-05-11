# Shows some PCE analysis from samples:
# - generation of samples
# - construction of PCE coefficients
# - sensitivity analysis of PCE coefficients
# - comparison of PCE predictor versus truth

import itertools
import numpy as np
from matplotlib import pyplot as plt

import indexing, wafp, sensitivity
from test_functions import genz_oscillatory
from opolynd import opolynd_eval
from recurrence import jacobi_recurrence

# construction of our "expensive model". Here just a collection of
# not-so-complicated test functions.
Nx = 100                    # Number of "spatial" grid points (here, space is 1D)
x = np.linspace(0, 1, Nx)   # spatial grid
print(x)

def expensive_model(zz,xx):
    """
    Evaluates a toy function as if it were an expensive model.

    x is a "spatial" grid (numpy array)
    z is a single d-variate parameter/random variable instance.

    Returns an array shaped like xx.
    """

    r = 2. # Fixed input into genz oscillatory function

    d = zz.size

    f = np.zeros(xx.shape)
    for ind,xval in np.ndenumerate(xx):
        f[ind] = genz_oscillatory(xval*np.ones(d), r)(zz)

    return f

######## Step 1: generate samples

d = 2                       # dimension of random/parameter space
k = 7                       # polynomial degree (parameter space)
poly_space = indexing.total_degree_indices

lambdas = poly_space(d, k)
N = lambdas.shape[0]

M = N + 10  # Number of samples (expensive model runs) to perform. Must
            # be at least N, but the +10 is arbitrary.
            # Of course, large M  ===> better.

candidate_mesh_size = int(2e2)

print("Generating parameter mesh...")
z = wafp.legendre_wafp(lambdas, M=candidate_mesh_size)
z = wafp.legendre_wafp_enrichment(z, lambdas, M-N)

# The samples are the array z, each row is a d-dimensional sample on the
# hypercube [-1,1]^d. The particular way these points are generated is
# random, so you'll get a different set of z each time this is fun, but
# the "quality" of the grid has relatively low variance from run to run.
# This takes a while, but these points may be stored for future use.

######## Step 2: run "expensive" model

u = np.zeros([Nx, M])

print("Evaluating model on mesh...")
for ind,zval in enumerate(z):
    u[:,ind] = expensive_model(zval, x)

print(u)

# Each column of u is a model run for a fixed parameter value

######## Step 3: compute PCE coefficients
print("Assembling PCE coefficients...")
ab = jacobi_recurrence(lambdas.max()+1, alpha=0., beta=0., probability=True)
V = opolynd_eval(z, lambdas, ab)
weights = np.sqrt(float(N)/float(M)) / np.sqrt(np.sum(V**2,axis=1))

# The PCE coefficients are computed as a weighted discrete least-squares
# estimator with the weights above.
coeffs = np.linalg.lstsq( (V.T*weights).T, (u*weights).T)[0].T

# Each row of coeffs contains PCE coefficients for a single x gridpoint.
# Each column of coeffs contains a particular PCE coefficient for all
# values of x.

######## Step 4: whatever postprocessing you want
print("Processing PCE coefficients...")
# Compute total sensitivities
total_sensitivities = sensitivity.pce_total_sensitivity(coeffs.T, lambdas, list(range(d)))

# Compute global sensitivities
# Compute main-effect and main-interaction sensitivities
Js = [[j] for j in range(d)]
[Js.append(comb) for comb in itertools.combinations(list(range(d)), 2)]

global_sensitivities = sensitivity.pce_global_sensitivity(coeffs.T, lambdas, Js)

# Compare surrogate discrepancy at validation points
M_validation = 100
z_validation = np.random.uniform(0, 1, [100, d])

u_truth = np.zeros([Nx, M_validation])
for ind, zval in enumerate(z_validation):
    u_truth[:,ind] = expensive_model(zval, x)

u_pce = np.zeros([Nx, M_validation])
V = opolynd_eval(z_validation, lambdas, ab)
u_pce = np.dot(V, coeffs.T).T

# Compute l2- and max-errors for each grid point:
l2_error = np.sqrt(np.sum((u_truth - u_pce)**2, axis=1))/np.sqrt(N)
linf_error = np.max(np.abs(u_truth - u_pce), axis=1)

print("L2 error on validation mesh: {0:1.3e}\nMaximum error on validation mesh: {1:1.3e}".format(np.linalg.norm(l2_error)/np.sqrt(Nx), np.max(linf_error.flatten())))

if d == 2:
    plt.plot(z[:,0], z[:,1], 'r.')
    plt.show()
