##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_label2annot
mris_label2annot_SOURCES=mris_label2annot.c
mris_label2annot_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_label2annot_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST = runit

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
