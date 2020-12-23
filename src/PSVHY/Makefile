########################################################################
#                                                                      #
# PSV Hybrid code makefile script                                      #
#                                                                      #
# Written by: Lianxing Wen, Seismological Lab.  Caltech.               #
# Last modified: Sat. 11, Aug. 1996.                                   #
#                                                                      #
########################################################################

BIN  = ../../bin/
GETPAR=../getpar
LIBS = -lm -L${GETPAR}/lib/ -lget
GETPARINC  = -I${GETPAR}/include
FC  = gfortran
CC   = gcc -Wall

all: psvfd aserpsvfd clean

genmodel: genmodel.o readmodel.o interpl.o nrutil.o messy_mod.o dsyssolve.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

psvfd:  psvfd.o readmodel.o interpl.o nrutil.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

psvfd.o: psvfd.c psvfd.h
	$(CC) -c psvfd.c $(GETPARINC)

aserpsvfd: aserpsvfd.o aser_psv.o readmodel.o aser.rout.o dcsource.o nrutil.o
	$(FC) -o $(BIN)/$@ $^ $(LIBS)

aserpsvfd.o: aserpsvfd.c
	$(CC) -c $(GETPARINC) $^ -o $@

clean:
	rm -f *.o
