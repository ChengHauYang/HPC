import numpy as np

def Function(x):
    return np.exp(-20*x**4)

def Function4(x):
    #return ((np.log(1-x))**2-x**2-x**3-11*x**4/12)/x**(5/2)
    return 5*x**(5/2)/6

def Function6(x):
    return np.sin(x**2)

def Function_test(x):
    return np.exp(3*x)*np.sin(2*x)

def CompTrap(f,a,b,n):
    h=(b-a)/n
    I=f(a)+f(b)
    for k in range(1,n):
        I=I+2*f(a+k*h)
    I=h*I/2
    return I

def CompSimps(f,a,b,n):
    h=(b-a)/n
    I=f(a)+f(b)+4*f(a+(n-1)*h)
    for k in range(1,int(n/2)):
        I=I+2*f(a+2*k*h)+4*f(a-h+2*k*h)
    I=h*I/3
    return I

def CompSimps38(f,a,b,n):
    h=(b-a)/n
    I=f(a)+f(b)+3*f(a+(n-1)*h)+3*f(a+(n-2)*h)
    for k in range(1,int(n/3)):
        I=I+2*f(a+3*k*h)+3*f(a-h+3*k*h)+3*f(a-2*h+3*k*h)
    I=h*3*I/8
    return I

def Romberg(f,a,b,n):
    # Romberg integration using n composite trapezoidal rule values
    h=b-a
    R=np.zeros((n,n))
    R[0][0]=h*(f(a)+f(b))/2
    for k in range(2,n+1):
        z=0
        for i in range(1,2**(k-2)+1):
            z=z+h*f(a+(i-1/2)*h)
        R[k-1][0]=(R[k-2][0]+z)/2
        h=h/2
        for j in range(2,k+1):
            R[k-1][j-1]=R[k-1][j-2]+(R[k-1][j-2]-R[k-2][j-2])/(4**(j-1)-1)
    return R

def print_R(R):
    n=len(R)
    for i in range(0,n):
        for j in range(0,i+1):
            print("%15.8e" %R[i][j],end=" ")
        print()

if __name__ == "__main__":

    a=-1
    b=1
    ExactSol=0.85722253700516515104
    Error=[]
    print("--------------composite Trapezoidal rule-------------")
    for k in range(1,12):
        n=2**(k+1)
        I=CompTrap(Function,a,b,n)
        Error.append(np.abs(I-ExactSol))
        print("k = %d, " %k,end="")
        print("n = %d, " %n,end="")
        print("approximations = %22.15e, " %I,end="")
        if (k>1):
            print("absolute errors = %22.15e, " %np.abs(I-ExactSol),end="")
            print("error ratios = %22.15e, " %(Error[k-2]/Error[k-1]),end="")
            print("order of convergence = %22.15e" %(np.log((Error[k-2]/Error[k-1]))/np.log(2)))
        else:
            print("absolute errors = %22.15e" %np.abs(I-ExactSol))


    Error=[]
    print("--------------composite Simpsons rule-------------")

    for k in range(1,12):
        n=2**(k+1)
        I=CompSimps(Function,a,b,n)
        Error.append(np.abs(I-ExactSol))
        print("k = %d, " %k,end="")
        print("n = %d, " %n,end="")
        print("approximations = %22.15e, " %I,end="")
        if (k>1):
            print("absolute errors = %22.15e, " %np.abs(I-ExactSol),end="")
            print("error ratios = %22.15e, " %(Error[k-2]/Error[k-1]),end="")
            print("order of convergence = %22.15e" %(np.log((Error[k-2]/Error[k-1]))/np.log(2)))
        else:
            print("absolute errors = %22.15e" %np.abs(I-ExactSol))


    Error=[]
    print("--------------composite Simpson’s three-eigths rule-------------")

    for k in range(1,12):
        n=3*2**k
        I=CompSimps38(Function,a,b,n)
        Error.append(np.abs(I-ExactSol))
        print("k = %d, " %k,end="")
        print("n = %d, " %n,end="")
        print("approximations = %22.15e, " %I,end="")
        if (k>1):
            print("absolute errors = %22.15e, " %np.abs(I-ExactSol),end="")
            print("error ratios = %22.15e, " %(Error[k-2]/Error[k-1]),end="")
            print("order of convergence = %22.15e" %(np.log((Error[k-2]/Error[k-1]))/np.log(2)))
        else:
            print("absolute errors = %22.15e" %np.abs(I-ExactSol))

    order=np.log(Error[10]/Error[0])/np.log(3*2/(3*2**11))
    print("order= %22.15e" %order)

    print("--------------composite Trapezoidal rule (P4 (c))-------------")
    n=30846
    a=0
    b=3/4
    I=CompTrap(Function4,a,b,n)
    print("approximations = %22.15e" %I)
    print("approximations (all) = %22.15e" %(I+1057/451))

    print("--------------Romberg-------------")
    n=0
    while 1:
        n+=1
        a=0
        b=2*np.sqrt(np.pi)
        ExactSol6=0.48624702331211065533
        R=Romberg(Function6,a,b,n)
        #print(R[n-1][n-1])
        if (abs(R[n-1][n-1]-ExactSol6))<1e-14:
            print(n)
            print_R(R)
            print()
            print_R(np.abs(R-ExactSol6))
            break

'''
    n=6
    a=0
    b=np.pi/4
    R=Romberg(Function_test,a,b,n)
    print(R[n-1][n-1])
    print_R(R)
'''
