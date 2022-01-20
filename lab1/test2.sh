
# bash script to check for software 
# Insert your name where indicated below, then type
#    $ bash test2.sh
# When it looks good, redirect output to a file:
#    $ bash test2h.sh > test2output.txt

echo 
echo Code run by  **insert your name here**
echo Environment variable ISUHPC is $ISUHPC
echo Environment variable MYHPC is $MYHPC

echo 
echo which ipython returns...
which ipython

echo 
echo which gcc returns...
which gcc

echo 
echo gcc --version returns...
gcc --version

echo
echo Compiling and running a C code...
gcc test3.c
./a.out

