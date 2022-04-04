CC = gcc
FFLAGS = -std=c99 -Wall -fopenmp
LFLAGS = -lgomp -lm

OBJECTS = main.o

main.exe:  $(OBJECTS)
	$(CC)  $(OBJECTS) -o main.exe $(LFLAGS)

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe
