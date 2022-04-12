CC = mpicc 
FFLAGS = -O3 -Wall -std=c99
LFLAGS = -lm
OBJECTS = main.o vector.o

main.exe:  $(OBJECTS)
	$(CC)  $(OBJECTS) -o main.exe $(LFLAGS)

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe *.data slurm*
