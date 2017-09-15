import numpy as np
from scipy import interpolate

from lowpassma import lowpassma

def tikhonov(A,L,Y,lambdas):
# Function to calculate the Greensite inverse solution
#   A - forward matrix
#   L - regularization matrix
#   Y - data
#   X_reg - inverse solution

    LL=(L.T@L)
    AA=(A.T@A)

    X_reg = np.zeros((AA.shape[0],Y.shape[1]))

    for i in range(0, Y.shape[1]):
    
        eta=np.zeros((lambdas.size,1))
        rho=np.zeros((lambdas.size,1))
        for k in range(0,lambdas.size):
            temp = np.linalg.solve(AA+lambdas[k]*LL, A.T@Y[:,i])
            eta[k]=np.linalg.norm(L*temp)
            rho[k]=np.linalg.norm(Y[:,i]-A*temp)         

        rholog=np.log10(rho)
        etalog=np.log10(eta)
        Trho=2*lambdas.size;

        tck = interpolate.splrep(rholog, etalog, s=0)
        etalog=interpolate.splev(np.linspace(np.min(rholog),np.max(rholog), Trho), tck, der=0);
        etalog=lowpassma(np.atleast_2d(etalog), 10);

        tck = interpolate.splrep(rholog, lambdas, s=0)
        lambdaspline=interpolate.splev(np.linspace(np.min(rholog),np.max(rholog), Trho), tck, der=0);

        detalog=np.diff(etalog, axis=1)
        signchanges=np.diff(np.sign(detalog))
        lastchangeind=np.where(signchanges<0)[-1]

        lamb = lambdaspline[lastchangeind]

        if not lamb:
            lamb = lambdas[int(np.ceil(lambdas.size/2))-1];

        X_reg[:,i] = np.linalg.solve(AA+lamb*LL,A.T@Y[:,i]);
    
        print('t={0}'.format(i))

    return X_reg
    
if __name__ == '__main__':
    print('Testing tikhonov')

    lambdas = 10**np.linspace(-2, 2, 1e3, endpoint=True)
    X = np.arange(0,100).reshape(5,20)    
    A = np.eye(5)
    Y = A@X
    L = np.eye(5)
    X_reg = tikhonov(A,L,Y,lambdas)

    print(X_reg)