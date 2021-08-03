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

extract_hjs.o: extract_hjs.f
	$(FOR) extract_hjs.f
extract_hjs: extract_hjs.o
	$(LINK) -o extract_hjs extract_hjs.o $(LIBSWIFT)
	rm extract_hjs.o
gen_hjs.o:	gen_hjs.f
	$(FOR) gen_hjs.f
gen_hjs: gen_hjs.o
	$(LINK) -o gen_hjs gen_hjs.o $(LIBSWIFT)
	rm gen_hjs.o
gen_hjs_tidal.o:	gen_hjs_tidal.f
	$(FOR) gen_hjs_tidal.f
gen_hjs_tidal: gen_hjs_tidal.o
	$(LINK) -o gen_hjs_tidal gen_hjs_tidal.o $(LIBSWIFT)
	rm gen_hjs_tidal.o
