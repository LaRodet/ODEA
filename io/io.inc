c***********************************************************************
c                        IO.INC
c***********************************************************************
c Include file for io routines

c...  Number of bytes in iflgchk
      integer IO_NBITS
      parameter(IO_NBITS=8)
c        bit 0 set ==> write big binary data file
c        bit 1 set ==> write real*8 binary file rather
c                       than real*4, ignored if bit0=F
c        bit 2 set ==> calc energy of system wrt time
c        bit 3 set ==> calc jacobi of the test particles
c        bit 4 set ==> check if particles are removed
c        bit 5 set ==> include J2 and J4 terms
c        bit 6 set ==> initial hierarchy check
c        bit 7 set ==> adaptive time step
c-----------------------------------------------------------------------
