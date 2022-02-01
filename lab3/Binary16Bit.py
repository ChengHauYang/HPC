def bin2int(bin_string):

    M=len(bin_string)-1

    # for checking (lab3)
    assert(isinstance(bin_string,str)) # check that bin string is a string
    assert(len(bin_string) == 16) # bin string has exactly 16 characters

    # check that each character is either ’0’ or ’1’
    a=[]
    for n in range(0,M+1):
        a.append(bin_string[n])
        if(int(a[n]) != 0 and int(a[n]) != 1):
            raise Exception('each character is either 0 or 1')

    a=[]
    x=0
    for n in range(0,M+1):
        a.append(bin_string[n])
        if n==0:
            x+=-int(a[n])*2**(M-n)
        else:
            x+=int(a[n])*2**(M-n)
    return x

def int2bin(x):

    # for checking (lab3)
    assert(type(x) is int) # check that n is an integer
    assert(x > -2**15-1) # not too small
    assert(x < 2**15) # not too large

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

def print_b2i(bin_string,i):
    print("Binary = %s, to be Integer %s"  %(bin_string,i))

def print_3b(bin1,bin2,binadd):
    print("Binary1 %s, + Binary2 %s = Added Binary %s"  %(bin1,bin2,binadd))

def print_i2b(bin_string,i):
    print("Integer = %s, to be Binary %s"  %(i,bin_string))

def print_i2b_wans(bin_string,hand_ans,i):
    print("Integer = %s, Calculated by hand = %s, Calculated by code %s"  %(i,hand_ans,bin_string))

def test():
    num=[0,-1,-32768,1,32767]
    ans=['0000000000000000','1111111111111111','1000000000000000','0000000000000001','0111111111111111']
    for i in range(len(num)):
        x = int2bin(num[i])
        print_i2b_wans(x,ans[i],num[i])
        if (x==ans[i]):
            print('-> coding ans matches hand ans')


def add2bins(bin1,bin2):

    M=len(bin1)-1

    a=[]
    b=[]
    x=0
    for n in range(0,M+1):
        a.append(bin1[n])
        b.append(bin2[n])

    sum=[]

    bin_out = []
    num_forward=0 # number that goes forward
    for n in range(0,M+1):
        index=M-n
        numer_now=num_forward+int(bin1[index])+int(bin2[index])

        if (numer_now ==3):
            num_forward=1
            sum.append(1)
        elif (numer_now ==2):
            num_forward=1
            sum.append(0)
        else:
            num_forward=0
            sum.append(numer_now)

    for n in range(0,M+1):
        index=M-n
        tmp = str(sum[index])
        bin_out.append(tmp)

    bin_out = ''.join(bin_out)
    return bin_out



if __name__ == "__main__":
    test()
    print("----------Add two Binary----------")    
    x=add2bins('0111111100000010','0111111100000010')
    print_3b('0111111100000010','0111111100000010',x)

'''
    print("----------Binary to Integer----------")
    x = bin2int('1111111100000010')
    print_b2i("1111111100000010",x)

    print("----------Integer to Binary----------")
    x = int2bin(-254)
    print_i2b(x,-254)
'''
