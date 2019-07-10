#! /bin/ksh
#
# Makefile to build simulations
#

FOR     =       gfortran -c -O3
LINK    =	gfortran

LIBGAG = -L$(GAG_EXEC_DIR)/lib -lgio -lgsys
SWIFT_DIR = .
ANAL = $(SWIFT_DIR)/anal
COORD = $(SWIFT_DIR)/coord
DISCARD = $(SWIFT_DIR)/discard
IO = $(SWIFT_DIR)/io
DRIFT = $(SWIFT_DIR)/drift
GETACCH = $(SWIFT_DIR)/getacch
KICKVH = $(SWIFT_DIR)/kickvh
STEP =  $(SWIFT_DIR)/step
ORBEL = $(SWIFT_DIR)/orbel
UTIL = $(SWIFT_DIR)/util
CE = $(SWIFT_DIR)/ce

clean:
	rm -f *.o *~ core

libswift:
	rm -f libswift.a
	$(FOR) $(CE)/*.f
	$(FOR) $(ANAL)/*.f
	$(FOR) $(COORD)/*.f
	$(FOR) $(DISCARD)/*.f
	$(FOR) $(IO)/*.f
	$(FOR) $(DRIFT)/*.f
	$(FOR) $(GETACCH)/*.f
	$(FOR) $(KICKVH)/*.f
	$(FOR) $(STEP)/*.f
	$(FOR) $(ORBEL)/*.f
	$(FOR) $(UTIL)/*.f
	ar -rv $(SWIFT_DIR)/libswift.a  *.o
	rm *.o
