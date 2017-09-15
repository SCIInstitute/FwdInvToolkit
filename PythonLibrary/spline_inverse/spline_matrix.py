import numpy as np
from scipy import interpolate
import scipy.sparse as sp

def splineMatrix(NumberOfKnotPoints,DimensionOfSpace,InterpolationDensity):
# function S = splineMatrix(NumberOfKnotPoints,DimensionOfSpace,InterpolationDensity)
# Author: Burak Erem
# 
# Given a matrix X of [DerivAtFirstKnot,KnotPoints,DerivAtLastKnot],
# this function produces a matrix S that multiplies a vectorized X:
# y = S*vec(X)
# Y=reshape(y,DimensionOfSpace,T) is a matrix of points on the interpolated
# spline curve where
# T=InterpolationDensity*(NumberOfKnotPoints-1)+NumberOfKnotPoints
# is the number of interpolated points (i.e. with InterpolationDensity
# number of points between each knot point)
# 
	N = NumberOfKnotPoints+2;
	S = sp.csr_matrix((InterpolationDensity*(NumberOfKnotPoints-1)+NumberOfKnotPoints,N), dtype=np.float32)

	for i in range(0, N):
		temp = np.zeros((N));
		temp[i] = 1;
		ss = interpolate.splrep(np.arange(1,N+1),temp); 
		tempcol = interpolate.splev(np.linspace(1,NumberOfKnotPoints,InterpolationDensity*(NumberOfKnotPoints-1)+NumberOfKnotPoints),ss);
		S[:,i] = tempcol.reshape((-1,1))

	S = sp.kron(S, sp.eye(DimensionOfSpace));

	return S

if __name__ == '__main__':
    print('testing spline_matrix')
    NumberOfKnotPoints = 100;
    DimensionOfSpace = 2;
    InterpolationDensity = 2;
    S = splineMatrix(NumberOfKnotPoints,DimensionOfSpace,InterpolationDensity)
    print(S)
