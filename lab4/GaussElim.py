def LUfactor(A):
    from numpy import array     # numpy array:   A = array([[....]])
    from numpy import shape     # compute the size of an array:  shape(A)
    from numpy import identity  # constructs the n x n identity matrix: identity(n)
    from numpy import copy
    #U=A # not correct
    U=copy(A)
    n=A.shape[0]
    L=identity(n)

    for k in range(1,n):
        for j in range(k+1,n+1):

            L[j-1][k-1]=U[j-1][k-1]/U[k-1][k-1]

            for i in range(k,n+1):
                U[j-1][i-1]=U[j-1][j-1]-L[j-1][k-1]*U[k-1][i-1]

            '''
            print(U[j-1][k-1])
            print(U[k-1][k-1])
            print(U[j-1][k-1]/U[k-1][k-1])
            '''
            
            
    return L,U
    
def determinant(A):
    [L,U]=LUfactor(A)
    
    detL=1
    detU=1
    n=A.shape[0]
    
    for i in range(n):
        detU=detU*U[i][i]
    
    detA=detU*detL
    
    return detA


if __name__ == "__main__":
    import numpy as np
    A1 = np.array([[5.0,-1.0,2.0],[10.0,3.0,7.0],[15.0,17.0,19.0]],dtype=float)
    [L1,U1]=LUfactor(A1)
    print('L1')
    print(L1)
    print('U1')
    print(U1)
    print('L1*U1')
    print(np.dot(L1,U1))
    
    print('--------------A1--------------')
    
    

    A2 = np.array([[4.0,1.0,0.0,0.0],[1.0,4.0,1.0,0.0],[0.0,1.0,4.0,1.0],[0.0,0.0,1.0,4.0]],dtype=float)
    [L2,U2]=LUfactor(A2)
    print('L2')
    print(L2)
    print('U2')
    print(U2)
    print('L2*U2')
    print(np.dot(L2,U2))

    print('--------------A2--------------')


    A3 = np.array([[1.0,-2.0,-2.0,-3.0],[3.0,-9.0,0.0,-9.0],[-1.0,2.0,4.0,7.0],[-3.0,-6.0,-26.0,2.0]],dtype=float)
    [L3,U3]=LUfactor(A3)
    print('L3')
    print(L3)
    print('U3')
    print(U3)

    print('--------------A3--------------')

    detA1=determinant(A1)
    #print(detA1)


