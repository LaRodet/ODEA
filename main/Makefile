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

# Simulation
odea.o: odea.f
	$(FOR) odea.f
odea: odea.o
	$(LINK) -fmax-stack-var-size=262144 -o odea odea.o  $(LIBSWIFT)
	rm odea.o
swift_hjs.o: swift_hjs.f
	$(FOR) swift_hjs.f
swift_hjs: swift_hjs.o
	$(LINK) -fmax-stack-var-size=262144 -o swift_hjs swift_hjs.o  $(LIBSWIFT)
	rm swift_hjs.o
odea_tidal.o: odea_tidal.f
	$(FOR) odea_tidal.f
odea_tidal: odea_tidal.o
	$(LINK) -fmax-stack-var-size=262144 -o odea_tidal odea_tidal.o  $(LIBSWIFT)
	rm odea_tidal.o
