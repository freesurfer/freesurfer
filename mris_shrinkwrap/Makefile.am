##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mris_shrinkwrap mris_AA_shrinkwrap
mris_AA_shrinkwrap_SOURCES=mris_AA_shrinkwrap.c
mris_AA_shrinkwrap_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_AA_shrinkwrap_LDFLAGS=$(OS_LDFLAGS)

mris_shrinkwrap_SOURCES=mris_shrinkwrap.c
mris_shrinkwrap_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_shrinkwrap_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
