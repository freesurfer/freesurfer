#
# Makefile for fsgdf
#
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

noinst_LIBRARIES= libfsgdf.a
libfsgdf_a_SOURCES = fsgdf.c

foodir=$(prefix)/lib
foo_DATA=libtclfsgdf$(DYNAMIC_LIB_EXT)

libtclfsgdf$(DYNAMIC_LIB_EXT): fsgdf_wrap.c fsgdf.c
	$(COMPILE) $(SHARED_CFLAGS) -c -o fsgdf_wrap.o fsgdf_wrap.c
	$(COMPILE) $(SHARED_CFLAGS) -DDECLARE_PROGNAME -c -o fsgdf.o fsgdf.c
	$(CC) $(SHARED_LDFLAGS) $(LDFLAGS) fsgdf.o fsgdf_wrap.o \
	$(tcl_libraries) \
	  $(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(LIBS) -o libtclfsgdf$(DYNAMIC_LIB_EXT)
	rm fsgdf.o

fsgdf_wrap.c: fsgdf.i
