CC = gcc
FC = gfortran
BIN  = ../../bin/

GETPAR = ../getpar
LIBS = -lm -L$(GETPAR)/lib -lget

all: aseriesCMB clean

OBJS = aseriesCMB.o aser_psv.o read_model.o aser.rout.o findcmbs.o \
	modify_model.o get_refls.o open_file.o diff.o

aseriesCMB: $(OBJS)
	$(FC) -o $(BIN)/$@ $^ $(LIBS)

aseriesCMB.o: aseriesCMB.c
	$(CC) -c $^ -o $@ -I$(GETPAR)/include

clean:
	rm -f *.o
