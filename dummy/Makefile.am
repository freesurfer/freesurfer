## 
## Makefile.am 
##
#######################################################################
# when you create your own program directory, copy the from to the to lines
# to your Makefile.am and change dummy to your program name
#
# 1. modify dev/configure.in to have your things
# 2. modify dev/Makefile.am to have your things
# 
# 3. automake (which will create Makefile.in)
# 4. autoconf (which will create configure)
# 5. configure .... and make to make sure that your new thing is working
# 6. checkin yourdir/Makefile.am, Makefile.in, yourprog.c or cpp.  (Don't checkin Makefile.)
# 7. checkin dev/configure.in, dev/configure, dev/Makefile.am, dev/Makefile.in
#
# if you have problems, please let me know (tosa@nmr.mgh.harvard.edu)
#########################################################################
# from ------------------------------------------------
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = dummy
dummy_SOURCES=dummy.c
dummy_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
dummy_LDFLAGS=$(OS_LDFLAGS)
# to --------------------------------------------------
#
# where $(OS_LDFLAGS) is to compile program statically for Linux
#
# by the way, you can have as many bin_PROGRAMS in the same directory
#
########################################################################
# when you want to create your program to test, use the following framework
# you can have as many programs (separated by space) as check_PROGRAMS
# your program can consist of many sources (you add to _SOURCES= line)
# You can do "make myown" to compile. The default CFLAGS="-g -O2", but
# you can do "make CFLAGS=-g myown" to create debugging version.
#
# the program listed as check_PROGRAMS are not compiled under standard "make".
# only "make check" will compile them.
# 
# please don't remove myown.c
########################################################################
check_PROGRAMS= myown

myown_SOURCES=myown.c
myown_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))



