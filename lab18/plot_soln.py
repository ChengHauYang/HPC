def plot_soln(create_png=False):

    import numpy as np; import string;
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    matplotlib.rcParams.update({'font.size': 16, 'font.family': 'sans-serif'});
    matplotlib.rcParams.update({'text.usetex': 'true'});

    # figure out how many processors were used
    fid = open('num_processors.data');
    linestring = fid.readline()
    linelist   = linestring.split()
    num_proc   = np.int64(linelist[0])
    fid.close()

    # figure out how many grid points there are
    fid = open('q0_0.data', 'r')
    mx = -1;
    while True:
        line = fid.readline()
        if not line: break
        mx = mx+1
    fid.close()

    x   = np.zeros((mx,num_proc,11),dtype=float)
    U   = np.zeros((mx,num_proc,11),dtype=float)
    V   = np.zeros((mx,num_proc,11),dtype=float)
    time = np.zeros(11,dtype=float);

    z = ['ro','g^','ro','g^','ro','g^','ro','g^','ro','g^','ro','g^','ro','g^','ro','g^'];

    for NumFrame in range(0,11):
        for NumProc in range(0,num_proc):

            # read-in all the floats
            fname = 'q' + str(NumFrame) + '_' + str(NumProc) + '.data'
            fid = open(fname, 'r')
            linestring = fid.readline()
            linelist   = linestring.split()
            time[NumFrame] = np.round(10000*np.float64(linelist[0]))/10000.0
            for k in range(0,mx):
                linestring    = fid.readline()
                linelist      = linestring.split()
                x[k,NumProc,NumFrame] = np.float64(linelist[0])
                U[k,NumProc,NumFrame] = np.float64(linelist[1])
                V[k,NumProc,NumFrame] = np.float64(linelist[2])
            fid.close()

    # exact solution
    Uex = np.zeros((mx*num_proc,11),dtype=float)
    Vex = np.zeros((mx*num_proc,11),dtype=float)
    xex = np.zeros((mx*num_proc,11),dtype=float)
    tmp_index=0
    for NumProc in range(0,num_proc):
        for k in range(0,mx):
            for NumFrame in range(0,11):
                xex[tmp_index,NumFrame] =  x[k,NumProc,NumFrame]
                Uex[tmp_index,NumFrame] =  func(x[k,NumProc,NumFrame],time[NumFrame])
            tmp_index = tmp_index+1

    # plot the result
    plt.rc("font", size=16);
    for NumFrame in range(0,11):

        plt.figure(NumFrame)

        title = "".join(("1D advection at time t = ",str(time[NumFrame])));
        plt.suptitle(title) # major title

        plt.subplot(1,2,1)
        plt.grid()
        plt.plot(xex[:,NumFrame],Uex[:,NumFrame],'b-',linewidth=1);

        for NumProc in range(0,num_proc):
            plt.plot(x[:,NumProc,NumFrame],U[:,NumProc,NumFrame],z[NumProc],markersize=4);

        plt.xlim(0.0,1.0)
        plt.xticks([0.0,0.25,0.5,0.75,1.0])

        plt.ylim(-0.1,1.1)
        plt.yticks([0.0,0.25,0.5,0.75,1.0])

        plt.subplot(1,2,2)
        plt.grid()
        plt.plot(xex[:,NumFrame],Uex[:,NumFrame],'b-',linewidth=1);

        for NumProc in range(0,num_proc):
            plt.plot(x[:,NumProc,NumFrame],V[:,NumProc,NumFrame],z[NumProc],markersize=4);


        plt.xlim(0.0,1.0)
        plt.xticks([0.0,0.25,0.5,0.75,1.0])

        plt.ylim(-0.1,1.1)
        plt.yticks([0.0,0.25,0.5,0.75,1.0])

        if(create_png):
            plt.savefig('advection'+str(NumFrame)+'.png', dpi=400, bbox_inches='tight')

    #plt.show()

def func(x,time):

    import numpy as np

    z = (x-time)-np.floor(x-time);
    f = np.exp(-400.0*pow(z-0.5,2))

    return f

if __name__ == "__main__":
    plot_soln(True);
