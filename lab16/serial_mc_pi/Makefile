CC = /usr/local/opt/llvm/bin/clang
FFLAGS = -I/usr/local/opt/llvm/include  -Wall -fopenmp
LFLAGS = -L/usr/local/opt/llvm/lib -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib  -lm
OBJECTS = main.o

main.exe:  $(OBJECTS)
	$(CC) $(LFLAGS) $(FFLAGS) $(OBJECTS) -o main.exe

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe
