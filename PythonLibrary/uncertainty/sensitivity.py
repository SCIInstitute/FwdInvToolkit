import numpy as np

def pce_total_sensitivity(coeffs, lambdas, js):
    # Computes total sensitivity associated to dimension js from PCE
    # coefficients. Coeffs may be an array, and the sensitivity is
    # computed across axis 0 (rows). js should be a list-type containing
    # dimension indices.
    # The output is len(js) x coeffs.shape[1]

    js = np.asarray(js,dtype=int)

    assert coeffs.shape[0] == lambdas.shape[0]
    assert np.all(js >= 0) and np.all(js < lambdas.shape[1])

    if len(np.shape(coeffs)) == 1:
        coeffs = np.reshape(coeffs, [coeffs.size, 1])

    unique_rows = np.vstack({tuple(row) for row in lambdas})
    # Just making sure
    assert unique_rows.shape[0] == lambdas.shape[0]

    variance_rows = np.linalg.norm(lambdas, axis=1) > 0.
    assert np.sum(np.invert(variance_rows)) == 1

    variances = np.sum(coeffs[variance_rows,:]**2, axis=0)
    total_sensitivities = np.zeros([js.size, coeffs.shape[1]])

    for (qj,j) in enumerate(js):
        total_sensitivities[qj,:] = np.sum(coeffs[lambdas[:,j]>0,:]**2, axis=0)/variances

    return total_sensitivities

def pce_global_sensitivity(coeffs, lambdas, Js):
    # Computes global sensitivity associated to dimensional indices js
    # from PCE coefficients. Coeffs may be an array, and the sensitivity is
    # computed across axis 0 (rows).
    # Js should be a list of index lists. The global sensitivity for each
    # index list is returned.
    # The output is len(Js) x coeffs.shape[1]

    assert coeffs.shape[0] == lambdas.shape[0]
    for j in Js:
        j = np.asarray(j)
        assert np.all(j >= 0) and np.all(j < lambdas.shape[1])

    if len(np.shape(coeffs)) == 1:
        coeffs = np.reshape(coeffs, [coeffs.size, 1])

    unique_rows = np.vstack({tuple(row) for row in lambdas})
    # Just making sure
    assert unique_rows.shape[0] == lambdas.shape[0]

    variance_rows = np.linalg.norm(lambdas, axis=1) > 0.
    assert np.sum(np.invert(variance_rows)) == 1

    variances = np.sum(coeffs[variance_rows,:]**2, axis=0)
    global_sensitivities = np.zeros([len(Js), coeffs.shape[1]])
    for (qj,j) in enumerate(Js):
        jc = np.setdiff1d(list(range(lambdas.shape[1])), j)
        inds = np.logical_and( np.all(lambdas[:,j] > 0, axis=1), \
                               np.all(lambdas[:,jc]==0, axis=1) )
        global_sensitivities[qj,:] = np.sum(coeffs[inds,:]**2, axis=0)/variances

    return global_sensitivities

if __name__ == "__main__":

    import indexing

    d = 2
    k = 2

    lambdas = indexing.total_degree_indices(d, k)

    M = 12

    coeffs = np.random.normal(size=[lambdas.shape[0], M])

    # Dimensional indices
    js = list(range(d))
    T = pce_total_sensitivity(coeffs, lambdas, js)

    # Sets of dimensional indices
    Js = [[0,], [1,], [0, 1]]
    G = pce_global_sensitivity(coeffs, lambdas, Js)
