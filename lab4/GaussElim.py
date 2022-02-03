def LUfactor(A):
    from numpy import array     # numpy array:   A = array([[....]])
    from numpy import shape     # compute the size of an array:  shape(A)
    from numpy import identity  # constructs the n x n identity matrix: identity(n)
    U=A
    n=A.shape[0]
    L=identity(n)

    for k in range(1,n):
        for j in range(k+1,n+1):
            L[j-1][k-1]=U[j-1][k-1]/U[k-1][k-1]
            for i in range(k,n+1):
                U[j-1][i-1]=U[j-1][j-1]-L[j-1][k-1]*U[k-1][i-1]
    return L,U


if __name__ == "__main__":
    import numpy as np
    A1 = np.array([[5,-1,2],[10,3,7],[15,17,19]])
    [L1,U1]=LUfactor(A1)
    print(L1)
    print(U1)
