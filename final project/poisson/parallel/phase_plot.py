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
    U = np.zeros(kmax,dtype=float)
    fid = open(string, 'r')
    num=0
    for k in range(0,kmax):
        linestring = fid.readline()
        linelist   = linestring.split()
        X[k]       = np.float64(linelist[0])
        Y[k]       = np.float64(linelist[1])
        U[k]       = np.float64(linelist[2])
        num+=1
    fid.close()

    # Third, return the result
    return X,Y,U,num;

def phase_plot():
    x,y,u,m2=getdata('output.data')

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.use('tkagg')
    from matplotlib import rc
    import numpy as np
    from matplotlib import pyplot, cm
    from mpl_toolkits.mplot3d import Axes3D
    import math
    
    rc('font',**{'family':'serif','serif':['Times']})
    rc('xtick', labelsize=18)
    rc('ytick', labelsize=18)

    fig,ax = plt.subplots(figsize=(8,6))
    plt.title("Solving Poisson Equation", fontsize=18)
    
    m=int(math.sqrt(m2))
    #print(u)
    u_OnGrid = np.zeros((m, m))
    x_OnGrid = np.zeros((m, m))
    y_OnGrid = np.zeros((m, m))
    #print(m)

    k=-1
    for i in range(0,m):
        for j in range(0,m):
            k+=1
            #print(i)
            #print(k)
            u_OnGrid[i,j]=u[k]
            x_OnGrid[i,j]=x[k]
            y_OnGrid[i,j]=y[k]
            
        


    pyplot.contourf(x_OnGrid, y_OnGrid, u_OnGrid, alpha=0.5, cmap='rainbow')
    pyplot.colorbar()
    #pyplot.contour(x_OnGrid, y_OnGrid, u_OnGrid)
    pyplot.xlabel('X')
    pyplot.ylabel('Y')
    pyplot.show()    
    ax.set_xlabel('u', fontsize=18)
    ax.set_ylabel('v', fontsize=18)
    plt.grid()
    #plt.show()


if __name__ == "__main__":
    phase_plot()
