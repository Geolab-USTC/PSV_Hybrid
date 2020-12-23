CC   = gcc -Wall
LIBS = -lm -L../getpar/lib -lget
FC  = gfortran
BIN  = ../../bin/
GETPARINC = -I../getpar/include

all: kirch Read_kirgreen clean

OBJS =	kirch.o readmodel.o convm.o convt.o conj.o  mins.o

kirch: $(OBJS)
	$(FC) -o $(BIN)/$@ $^ $(LIBS)

Read_kirgreen: Read_kirgreen.o diff.o readmodel.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

kirch.o: kirch.c
	$(CC) -c kirch.c $(GETPARINC)

Read_kirgreen.o: Read_kirgreen.c
	$(CC) -c Read_kirgreen.c ${GETPARINC}

clean:
	rm -f *.o
