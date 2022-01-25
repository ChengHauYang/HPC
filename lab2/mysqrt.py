''' my own square root function '''

def sqrt_newton(x,s=1.0,printshow=False):

    kmax = 10;

    for k in range(0,kmax):
        if printshow:
            print("Before iteration %s, s = %20.15f" % (k,s));
        s = 0.5*(s + x/s);

    if printshow:
        print("After %s iterations, s = %20.15f" % (k+1,s));

    return s;

def sqrt_secant(x,s=1.0,s_minus=0.0,printshow=False):

    kmax = 10;

    for k in range(0,kmax):
        if printshow:
            print("Before iteration %s, s = %20.15f" % (k,s));
        s = (x+s*s_minus)/(s+s_minus);
        s_minus=s;

    if printshow:
        print("After %s iterations, s = %20.15f" % (k+1,s));

    return s;

if __name__ == "__main__":

    print("----------Newton method----------")
    s = sqrt_newton(15.0,1.0,True);

    print("----------Secant method----------")
    s = sqrt_secant(15.0,1.0,0.0,True);

#    print(" ")

#    s = sqrt(8.0,7.0,True);
