import numpy as np

def PairwiseDistance(input, maxcolsize=[]):
# PairwiseDistance This function calculates pairwise distances (2norm SQUARED) between points
# in X and Y or, if this function is called with one argument, between
# points in X only.
# Author: Ramon Martinez Orellana
# Updated for large inputs: Burak Erem
# 
# D = PairwiseDistance(A,[B])
#   INPUT ARGS
#       A: matrix with points as columns
#       B: (optional) matrix with points as columns
#   OUTPUT ARGS (NOTE: "distance" here means 2-norm SQUARED of the difference between points)
#       D: distance matrix. D(i,j) corresponts to the distance between the
#       ith point in A and the jth point in B. If
#       there is only one argument, D(i,j) is the distance between ith
#       point in A and jth point in B.
# 
    nargin = len(input)

    if not maxcolsize:
        if nargin == 1:
            # Calculates the pairwise distances squared between columns of a matrix X
            # dist(xi,xj) = (xi - xj)^2 = ||xi||^2 + ||xj||^2 - 2*xi'*xj
            # //2norm squared
            D = distSingle(input[0]);
    
        if nargin == 2:
            # Calculates the pairwise distances between points in matrix X and Y
            # dist(xi,yj) = (xi - yj)^2 = ||xi||^2 + ||yj||^2 - 2*xi'*yj
            # 2norm squared
            D = distDouble(input[0],input[1]);
        
    else:
        # Block-column-wise computations        
        N1=input[0].shape[1]
        N2=input[1].shape[1]
        
        blocks1=np.arange(0,N1,maxcolsize)
        if np.max(blocks1) != N1:
            blocks1 = np.append(blocks1, N1)

        blockcount1=blocks1.size
        
        blocks2=np.arange(0,N2,maxcolsize);
        if np.max(blocks2) != N2:
            blocks2 = np.append(blocks2, N2)

        blockcount2=blocks2.size;
        
        D = np.zeros((N1,N2));
        
        for i in range(0,blockcount1-1):
            blockrange1=np.arange(blocks1[i],blocks1[i+1]);
            for j in range(0,blockcount2-1):
                blockrange2=np.arange(blocks2[j],blocks2[j+1])
                
                D[np.ix_(blockrange1,blockrange2)] = distDouble(input[0][:,blockrange1],input[1][:,blockrange2])    

    return D

def distSingle(singleMatrix):
	K = singleMatrix.T@singleMatrix;
	vD = np.diag(K);
	N = vD.size;
	A = np.tile(vD.reshape(-1,1),(1,N))
	D = -2*K + A + A.T
	return D

def distDouble(firstMatrix,secondMatrix):
	Dy = np.tile(np.sum(secondMatrix**2,axis=0),(firstMatrix.shape[1],1))
	Dx = np.tile(np.sum(firstMatrix**2,axis=0).reshape(-1,1),(1,secondMatrix.shape[1]))
	Dxy = firstMatrix.T@secondMatrix
	D = Dy + Dx - 2*Dxy
	return D


if __name__ == '__main__':

    A = np.arange(0,100).reshape(20,5)
    D = PairwiseDistance([A])
    print(D)
    B = np.arange(0,100).reshape(20,5)+50
    D = PairwiseDistance([A, B])
    print(D)
    D = PairwiseDistance([A, B], 3)
    print(D)