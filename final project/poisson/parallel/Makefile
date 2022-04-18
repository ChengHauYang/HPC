# $ISUHPC/lectures/lecture12/codes/matrix/Makefile
CC = mpicc
FFLAGS = -O3 -Wall
LFLAGS = -lm
OBJECTS = main.o \
matrix.o 

.PHONY: clean

main.exe: $(OBJECTS)
	$(CC)  $(OBJECTS) -o main.exe $(LFLAGS)

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) main.exe
