#CC = /usr/local/opt/llvm/bin/clang
#FFLAGS = -I/usr/local/opt/llvm/include -Wall -fopenmp
#FLAGS = -L/usr/local/opt/llvm/lib -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
CC = /usr/local/opt/llvm/bin/clang
FFLAGS = -I/usr/local/opt/llvm/include -std=c99  -Wall 
LFLAGS = -L/usr/local/opt/llvm/lib -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib  -lm

OBJECTS = main.o

main.exe:  $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main.exe

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe
