def read_input_file():
    '''
    This function reads the file "input.data"
    and stores the results in the NumPy array A.
    '''

    # First, figure out how many floats there are
    fid = open('input.data', 'r')
    kmax = 0;
    while True:
        line = fid.readline()
        if not line: break
        kmax = kmax+1
    fid.close()

    import numpy as np

    # Second, read-in all the floats
    A = np.zeros(kmax,dtype=float)
    fid = open('input.data', 'r')
    for k in range(0,kmax):
        linestring = fid.readline()
        linelist   = linestring.split()
        A[k]       = np.float64(linelist[0])
    fid.close()

    # Third, return the result
    return A;

if __name__ == "__main__":
    A = read_input_file()
    print ("Number of elements read = ",A.shape[0])
