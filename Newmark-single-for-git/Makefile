#!/bin/csh
#
CMP = gfortran
BIN = /Vvn_Codes/bin 
SRC = ./src/
SRCP = ./src/
#
MAIN = $(SRC)Newmark-single.f
OBJ = Newmark-single.o  BPFILTER.o 
#
#
Newmark-single.o: $(MAIN)
	$(CMP) $(MAIN)  -c

BPFILTER.o: $(SRC)BPFILTER.f
	$(CMP) $(SRC)BPFILTER.f -c

Newmark-single: $(OBJ)
	$(CMP) $(OBJ) -o $@
#	cp Newmark-single $(BIN)
#

clean:
	rm -f *.o *~ *.ps

