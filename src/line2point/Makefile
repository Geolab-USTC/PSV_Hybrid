FC  = gfortran -Wall
CC   = gcc -Wall

CFLAGS = -I../getpar/include

BIN  = ../../bin/

LIBS = -L$(SACHOME)/lib -lsac -lsacio -L../getpar/lib -lm  -lget

all: line2point clean

line2point: readhead.o line2point.o diff.o
	$(FC) -o $(BIN)/$@ $^ $(LIBS)

clean:
	rm -f *.o
