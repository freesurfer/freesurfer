#
# Makefile.am for libutils.a
#
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

lib_LIBRARIES= libutils.a
libutils_a_SOURCES = $(wildcard *.c)
am_libutils_a_OBJECTS = $(libutils_a_SOURCES:%.c=%.o)

# tqli needs with -ffloat-store
tqli.o: tqli.c
	@echo "compiling tqli.c to tqli.o (not optimized)"
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
	$(CPPFLAGS) $(AM_CFLAGS) -ffloat-store -c `test -f '$<' || echo '$(srcdir)/'`$<

.c.o:
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
	$(CPPFLAGS) $(AM_CFLAGS) -c `test -f '$<' || echo '$(srcdir)/'`$<

