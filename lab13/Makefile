# $ISUHPC/lectures/lecture12/codes/matrix/Makefile
CC = gcc
FFLAGS = -O3 -Wall
LFLAGS = 
OBJECTS = main.o \
matrix.o \
trimatrix.o \
Hessenberg.o \
QRA.o \
QRdecomp.o \
FindZero.o \
GetShift.o

main.exe: $(OBJECTS)
	$(CC)  $(OBJECTS) -o main.exe $(LFLAGS)

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) main.exe
