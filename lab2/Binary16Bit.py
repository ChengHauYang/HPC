def bin2int(bin_string):
    a=[]
    M=len(bin_string)-1
    x=0
    for n in range(0,M+1):
        a.append(bin_string[n])
        if n==0:
            x+=-int(a[n])*2**(M-n)
        else:
            x+=int(a[n])*2**(M-n)
    return x

def int2bin(x):
    bin_string = [];
    if x==0:
        for n in range(0,16):
            bin_string.append('0')
    else:
        if x<0:
            tmp = '1'
            xmod = x + 2**15
        else:
            tmp = '0'
            xmod = x
        bin_string.append(tmp)

        for n in range(1,16):
            pw = 2**(15-n)
            z = int(xmod/pw) # fix for python3: int/int=float(py3)
            xmod = xmod - z*pw
            tmp = str(z)
            bin_string.append(tmp)
    bin_string = ''.join(bin_string)
    return bin_string

if __name__ == "__main__":

    print("----------Binary to Integer----------")
    x = bin2int('1111111100000010')
    print("1111111100000010")
    print("to")
    print(x)

    print("----------Integer to Binary----------")
    x = int2bin(-254)
    print(-254)
    print("to")
    print(x)
