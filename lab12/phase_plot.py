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

def phase_plot():
    x,y=getdata('cputime.data')

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.use('tkagg')
    from matplotlib import rc
    import numpy as np

    rc('font',**{'family':'serif','serif':['Times']})
    rc('xtick', labelsize=18)
    rc('ytick', labelsize=18)

    fig,ax = plt.subplots(figsize=(8,6))
    plt.title("CPU Time Versus Matrix Size", fontsize=18)
    #ax.plot(x,y,color='red')
    plt.semilogy(x, y, marker="o", linewidth=2)


    ax.set_xlabel('Matrix Size NxN', fontsize=18)
    ax.set_ylabel('CPU Time (sec)', fontsize=18)
    plt.grid()
    #plt.show()
    plt.savefig('cputime.png', bbox_inches='tight')


if __name__ == "__main__":
    phase_plot()
