def getdata(string):
    '''
         INPUT:
         ------
         infile is a string and refers to a file
         that has 2 columns, the first column
         has the x-coordinates and the second
         column has the y-coordinates of the
         input data.
         WHAT THIS FUNCTION DOES:
         ------------------------
         This function reads in all the data
         from the file "infile" and returns the
         data as 2 vectors (numpy arrays): x and y.
         OUPUT:
         ------
         2 numpy arrays (vectors): x and y.
    '''

    # First, figure out how many floats there are
    fid = open(string, 'r')
    kmax = 0;
    while True:
        line = fid.readline()
        if not line: break
        kmax = kmax+1
    fid.close()

    import numpy as np

    # Second, read-in all the floats
    X = np.zeros(kmax,dtype=float)
    Y = np.zeros(kmax,dtype=float)
    fid = open(string, 'r')
    for k in range(0,kmax):
        linestring = fid.readline()
        linelist   = linestring.split()
        X[k]       = np.float64(linelist[0])
        Y[k]       = np.float64(linelist[1])
    fid.close()

    # Third, return the result
    return X,Y;

def approx(x,y):
    '''
      INPUT:
      ------
      x and y are 2 numpy arrays (vectors) of
      the same length.
      WHAT THIS FUNCTION DOES:
      ------------------------
      This function takes the data (x,y) and calculates
      the 3 coefficients, a[0], a[1], and a[2], of a quadratic function:
             y = a[0] + a[1]*x + a[2]*x**2
      The coefficients are computed by forcing the quadratic
      polynomial to give the best-fit to the input data in
      the least-squares sense.
      OUTPUT:
      -------
      numpy array (vector) of length 3: a.
    '''
    import numpy as np

    # set up A
    #print(len(x))
    A=np.zeros([len(x), 3])
    B=np.zeros([len(x), 1])
    for i in range(len(x)):
        A[i][0]=1
        A[i][1]=x[i]
        A[i][2]=x[i]**2
        B[i]=y[i]
    #print(A)
    #print(B)
    #a=np.linalg.solve(A, B)
    #print(a)
    #print('=============')

    #normal equations
    RHS=np.dot(A.transpose(),A)
    LHS=np.dot(A.transpose(),B)
    #print(RHS)
    #print(LHS)
    a=np.linalg.solve(RHS, LHS)
    #print(np.allclose(np.dot(RHS, a), LHS)) # checking

    return a

def plotdata(x,y,a,infile):
    '''
      INPUT:
      ------
      x and y are 2 numpy arrays (vectors) of
      the same length. a is a numpy array (vector)
      of length 3. infile is a string that is the
      name of the input data file.
      WHAT THIS FUNCTION DOES:
      ------------------------
      This function takes the data (x,y) and
      the 3 coefficients, a[0], a[1], and a[2], of a quadratic function:
             y = a[0] + a[1]*x + a[2]*x**2
      and plots them on the same graph.
      OUTPUT:
      -------
      None. A PNG image file is created in the current directory.
    '''
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.use('tkagg')
    from matplotlib import rc
    import numpy as np

    rc('font',**{'family':'serif','serif':['Times']})
    rc('xtick', labelsize=18)
    rc('ytick', labelsize=18)

    fig,ax = plt.subplots(figsize=(8,6))
    plt.title(infile + ": Data + Least Squares Fit w/ Quad. Poly.", fontsize=18)
    xx=np.linspace(min(x), max(x), num=100)
    ax.plot(xx,a[0] + a[1]*xx + a[2]*xx**2,'--',color='red')
    plt.scatter(x,y)
    plt.legend(['quadratic approximation','data points'], prop={"size":15})

    ax.set_xlabel('x', fontsize=18)
    ax.set_ylabel('y', fontsize=18)
    plt.grid()
    #plt.show()
    plt.savefig(infile+'.png', bbox_inches='tight')


if __name__ == "__main__":
    X,Y=getdata('xydata1.dat')
    print("x=%s" %X)
    print("y=%s" %Y)
    a=approx(X,Y)
    print("a=%s" %a)
    plotdata(X,Y,a,'xydata1.dat')
    print('======================')

    X,Y=getdata('xydata2.dat')
    print("x=%s" %X)
    print("y=%s" %Y)
    a=approx(X,Y)
    print("a=%s" %a)
    plotdata(X,Y,a,'xydata2.dat')
    print('======================')

    X,Y=getdata('xydata3.dat')
    print("x=%s" %X)
    print("y=%s" %Y)
    a=approx(X,Y)
    print("a=%s" %a)
    plotdata(X,Y,a,'xydata3.dat')
    print('======================')

    X,Y=getdata('xydata4.dat')
    print("x=%s" %X)
    print("y=%s" %Y)
    a=approx(X,Y)
    print("a=%s" %a)
    plotdata(X,Y,a,'xydata4.dat')
    print('======================')
