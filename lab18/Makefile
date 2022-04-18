# $ISUHPC/lectures/lecture23/codes/MPI_Advection1D/Makefile
CC = mpicc
FFLAGS = -O3 -Wall 
LFLAGS = -lm
OBJECTS = vector.o main.o

main.exe:  $(OBJECTS)
	$(CC) $(OBJECTS) -o main.exe  $(LFLAGS)

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe *.png *.data slurm*
