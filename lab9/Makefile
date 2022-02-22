# $ISUHPC/lectures/lecture11/codes/fraction/Makefile6
CC = gcc
FFLAGS = -O3
LFFLAG =
OBJECTS = main.o \
	fraction_reduce.o get_prime_factors.o\
	fraction_add.o fraction_substract.o fraction_mult.o fraction_divide.o\
	make_proper.o

.PHONY: clean help

output.txt: main.exe
	./main.exe > output.txt

main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main.exe

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) main.exe output.txt

help:
	@echo "Valid targets:"
	@echo "  main.exe"
	@echo "  main.o"
	@echo "  fraction_add.o"
	@echo "  fraction_substract.o"
	@echo "  fraction_mult.o"
	@echo "  fraction_divide.o"
	@echo "  fraction_reduce.o"
	@echo "  get_prime_factors.o"
	@echo "  make_proper.o"
	@echo "  clean:  removes *.o, *.txt, and *.exe files"
