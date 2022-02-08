def quicksort(Ain,lo=0,hi=-100):
    '''
This function uses the quicksort algorithm
to sort the input NumPy array Ain and return
a NumPy array Aout, which is the sorted version
of Ain (sorted from lowest value to greatest)
    '''
    lo=lo-1
    hi=hi-1
    if (hi==-100):
        hi = len(Ain)-1

    if lo < hi:
        p = partition(Ain, lo, hi)
        quicksort(Ain, lo + 1, p)
        quicksort(Ain, p + 2, hi + 1)

    Aout=Ain
    return Aout

def partition(A, lo, hi):
    pivot = A[hi]
    i = lo # place for swapping
    for j in range(lo,hi): # lo to hi - 1
        if (A[j] <= pivot): # swap A[i] with A[j]
            temp=A[i]
            A[i]=A[j]
            A[j]=temp
            i=i+1
    # swap A[i] with A[hi]
    temp2=A[i]
    A[i]=A[hi]
    A[hi]=temp2
    return i


def filesort():
    '''
      This function does 3 things:
          1. reads the file input.data, which contains a random
              number of of random floats, into a NumPy array.
          2. Uses quicksort to sort the input.data file
          3. Output the sorted data into a new file called sorted.data.
    '''

    A = read_input_file()
    print(A)
    Aout=quicksort(A, 1, len(A))
    generate_output_file(Aout)
    return Aout


def generate_output_file(A):

    fid = open('input.data', 'w')
    for k in range(0,len(A)):
        value = rand.uniform(0.0, 1000.0)
        fid.write("%12.6e" % A[k])
        fid.write("\n");
    fid.close()

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
    import numpy as np
    A1 = np.array([1.0, 2.0, 3.0])
    Asort1=quicksort(A1, 1, len(A1))
    print(Asort1)

    A2 = np.array([3.0, 2.0, 1.0])
    Asort2=quicksort(A2, 1, len(A2))
    print(Asort2)

    A3 = np.array([-3.0, 2.0, -1.0])
    Asort3=quicksort(A3, 1, len(A3))
    print(Asort3)

    A4 = np.array([-3.0, 2.0, -1.0])
    Asort4=quicksort(A4, 1, len(A4))
    print(Asort4)

    A5 = np.array([1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0])
    Asort5=quicksort(A5, 1, len(A5))
    print(Asort5)


    filesort()
'''
    A6 = np.array([10, 7, 8, 9, 1, 5])
    Asort6=quicksort(A6, 1, len(A6))
    print(Asort6)
'''
