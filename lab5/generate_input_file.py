def generate_input_file():
    '''
    This function creates a file called
    input.data. This file will have a
    random number of floats each of which
    is a random number between 0.0 and 1000.0
    '''
    
    import random as rand
    kmax = rand.randint(1,1000)

    fid = open('input.data', 'w')    
    for k in range(0,kmax):
        value = rand.uniform(0.0, 1000.0)
        fid.write("%12.6e" % value)
        fid.write("\n");
    fid.close()

if __name__ == "__main__":
    generate_input_file()
