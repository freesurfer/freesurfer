## 
## Makefile.am 
##

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_parselabel
mri_parselabel_SOURCES=mri_parselabel.cpp
mri_parselabel_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_parselabel_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
