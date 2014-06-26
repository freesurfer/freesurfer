##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_parcellate_connectivity
mris_parcellate_connectivity_SOURCES=mris_parcellate_connectivity.c
mris_parcellate_connectivity_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_parcellate_connectivity_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
