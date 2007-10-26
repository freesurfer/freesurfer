#
# Makefile.am
#

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = label_border
label_border_SOURCES=label_border.c
label_border_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
label_border_LDFLAGS=$(OS_LDFLAGS)

#check_PROGRAMS= myown
#myown_SOURCES=myown.c
#myown_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
