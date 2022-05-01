# $ISUHPC/lectures/lecture12/codes/matrix/Makefile
#CC = mpicxx
CXX = mpicxx
FFLAGS = -O3 -Wall
LFLAGS = -lm
OBJECTS = main.o \
matrix.o 

.PHONY: clean

main.exe: $(OBJECTS)
	$(CXX)  $(OBJECTS) -o main.exe $(LFLAGS)

%.o: %.c
	$(CXX) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) main.exe
