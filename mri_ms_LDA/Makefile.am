##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_ms_LDA
mri_ms_LDA_SOURCES=mri_ms_LDA.c
mri_ms_LDA_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_ms_LDA_LDFLAGS=$(OS_LDFLAGS)

#check_PROGRAMS=myown
#myown_SOURCES=myown.c
#myown_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
