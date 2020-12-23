CC   = gcc -Wall
CFLAGS = -I../getpar/include
LIBS = -lm -L../getpar/lib/ -lget
BIN  = ../../bin

all:  demult demult2kir clean

demult: demult.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

demult2kir: demult2kir.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

clean:
	rm -f *.o
