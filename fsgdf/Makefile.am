#
# Makefile for fsgdf
#
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

noinst_LIBRARIES= libfsgdf.a
libfsgdf_a_SOURCES = fsgdf.c

# until MacOSX is fixed we don't build shared lib
# foodir=$(prefix)/lib
# foo_DATA=libtclfsgdf$(DYNAMIC_LIB_EXT)
#
# MacOSX requires to tell the libs
libtclfsgdf$(DYNAMIC_LIB_EXT): fsgdf_wrap.c fsgdf.c
	$(COMPILE) $(SHARED_CFLAGS) -c -o fsgdf_wrap.o fsgdf_wrap.c
	$(COMPILE) $(SHARED_CFLAGS) -DDECLARE_PROGNAME -c -o fsgdf.o fsgdf.c
	$(CC) $(SHARED_LDFLAGS) $(LDFLAGS) fsgdf.o fsgdf_wrap.o \
	-o libtclfsgdf$(DYNAMIC_LIB_EXT) \
	rm fsgdf.o

fsgdf_wrap.c: fsgdf.i
