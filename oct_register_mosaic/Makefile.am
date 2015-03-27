##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = oct_register_mosaic 
oct_register_mosaic_SOURCES=oct_register_mosaic.c
oct_register_mosaic_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
oct_register_mosaic_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
