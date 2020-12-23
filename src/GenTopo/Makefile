CFLAGS = -Wall
BIN  = ../../bin

all: gentopo clean

gentopo: gentopo.o
	$(CC) $^ -o ${BIN}/$@ -lm

clean:
	rm *.o
