import numpy as np

def genz_oscillatory(c, r):
    """
    Returns a multidimensional function evaluating the Genz
    `oscillatory' test function, defined as

        f(x) = cos( 2*pi*r + sum( dot(x, c) ) )

    where x has the same column dimension as c.

    c is a vector, and r is a scalar. The dimension of the input to the
    function is inferred from the size of c.
    """

    return lambda x: np.cos( 2*np.pi*r + np.dot(x, c) )

def genz_oscillatory_integral(c, r):
    """
    When c is a vector of size d, returns the value of the integral

      \int__{[-1,1]^d} f(x) 1/2^d dx,

    where f(x) is the function defined in genz_oscillatory.
    """

    return np.prod(np.sinc(c/np.pi)) * np.cos(2*np.pi*r)
