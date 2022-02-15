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
    u1,v1=getdata('output1.data')
    u2,v2=getdata('output2.data')
    u3,v3=getdata('output3.data')
    u4,v4=getdata('output4.data')
    u5,v5=getdata('output5.data')

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.use('tkagg')
    from matplotlib import rc
    import numpy as np

    rc('font',**{'family':'serif','serif':['Times']})
    rc('xtick', labelsize=18)
    rc('ytick', labelsize=18)

    fig,ax = plt.subplots(figsize=(8,6))
    plt.title("Solving ODE using Verlet Scheme", fontsize=18)
    ax.plot(u1,v1,color='red')
    ax.plot(u2,v2)
    ax.plot(u3,v3)
    ax.plot(u4,v4)
    ax.plot(u5,v5)
    plt.legend(['IC 1','IC 2','IC 3','IC 4','IC 5'], prop={"size":15})

    ax.set_xlabel('u', fontsize=18)
    ax.set_ylabel('v', fontsize=18)
    plt.grid()
    #plt.show()
    plt.savefig('verlet.png', bbox_inches='tight')


if __name__ == "__main__":
    phase_plot()
