CC = clang -Wall
BIN = ../../bin
GETPAR = ../getpar
GETPARLIB = $(GETPAR)/lib
GETPARINC = $(GETPAR)/include
CFLAGS = -I$(GETPARINC)

all: check_wkm check_grt check_fdkir check_grtkir check_kirconv \
	check_fdtrace

check_wkm: check_wkm.o sacio.o
	$(CC) -o $(BIN)/$@ $^ -L$(GETPARLIB) -lget

check_grt: check_grt.o sacio.o
	$(CC) -o $(BIN)/$@ $^ -L$(GETPARLIB) -lget

check_fdkir: check_fdkir.o sacio.o
	$(CC) -o $(BIN)/$@ $^ -L$(GETPARLIB) -lget

check_grtkir: check_grtkir.o sacio.o
	$(CC) -o $(BIN)/$@ $^ -L$(GETPARLIB) -lget

check_kirconv: check_kirconv.o sacio.o
	$(CC) -o $(BIN)/$@ $^ -L$(GETPARLIB) -lget

check_fdtrace: check_fdtrace.o sacio.o
	$(CC) -o $(BIN)/$@ $^ -L$(GETPARLIB) -lget

clean:
	rm *.o

