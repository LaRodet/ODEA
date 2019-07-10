#! /bin/ksh
#
# Makefile to build simulations
#

FOR     =       gfortran -c -O3
LINK    = gfortran
SWIFT_DIR = ..
LIBSWIFT  = -L$(SWIFT_DIR) -lswift

clean:
	rm -f *.o *~ core

extract_hjs.o:	extract_hjs.f
	$(FOR) extract_hjs.f
extract_hjs: extract_hjs.o
	$(LINK) -o extract_hjs extract_hjs.o $(LIBSWIFT)
gen_hjs.o:	gen_hjs.f
	$(FOR) gen_hjs.f
gen_hjs: gen_hjs.o
	$(LINK) -o gen_hjs gen_hjs.o $(LIBSWIFT)
