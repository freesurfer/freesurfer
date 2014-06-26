##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mris_deform mris_ca_deform mris_nudge
mris_deform_SOURCES=mris_deform.c
mris_deform_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_deform_LDFLAGS=$(OS_LDFLAGS)

mris_ca_deform_SOURCES=mris_ca_deform.c
mris_ca_deform_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_ca_deform_LDFLAGS=$(OS_LDFLAGS)

mris_nudge_SOURCES=mris_nudge.c
mris_nudge_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_nudge_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

