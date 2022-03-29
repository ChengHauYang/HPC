CC = /usr/local/opt/llvm/bin/clang
FFLAGS = -I/usr/local/opt/llvm/include -Wall -fopenmp
LFLAGS = -L/usr/local/opt/llvm/lib -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
OBJECTS = main_local.o matrix_local.o

main.exe:  $(OBJECTS)
	$(CC) $(LFLAGS) $(FFLAGS) $(OBJECTS) -o main_local.exe

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe
