#
# Makefile.am for NrrdIO directory
#

AM_CPPFLAGS=-I$(top_srcdir)/include/NrrdIO \
	-DTEEM_DIO=$(TEEM_DIO) -DTEEM_32BIT=$(TEEM_32BIT) \
	-DTEEM_QNANHIBIT=$(TEEM_QNANHIBIT) -DTEEM_ENDIAN=$(TEEM_ENDIAN)

noinst_LIBRARIES= libnrrdio.a

libnrrdio_a_SOURCES = $(wildcard *.c) $(wildcard *.cpp)

am_libnrrdio_a_OBJECTS = $(libnrrdio_a_SOURCES:%.c=%.o) \
	$(libnrrdio_a_SOURCES:%.cpp=%.o)

.c.o:
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

.cpp.o:
	@echo "compiler: $(COMPILE) end_compiler"
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
