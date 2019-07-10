#! /bin/ksh
#
# Makefile to build simulations
#

FOR     =       gfortran -c -O3 -ffpe-trap=invalid,zero,overflow,underflow
CC      =       cc -c
LINK = gfortran
SWIFT_DIR = ..
LIBSWIFT  = -L$(SWIFT_DIR) -lswift

clean:
	rm -f *.o *~ core

special.o:	special.f
	xlf -c special.f

# Simulation
odea.o: odea.f
	$(FOR) odea.f
odea: odea.o
	$(LINK) -fmax-stack-var-size=262144 -o odea odea.o  $(LIBSWIFT)
	rm odea.o
