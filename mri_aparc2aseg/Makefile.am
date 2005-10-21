## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_aparc2aseg
mri_aparc2aseg_SOURCES=mri_aparc2aseg.c
mri_aparc2aseg_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_aparc2aseg_LDFLAGS=$(OS_LDFLAGS)


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
