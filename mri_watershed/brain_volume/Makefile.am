## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_brain_volume
mri_brain_volume_SOURCES=mri_brain_volume.cpp mrivariables.h TVector.h
mri_brain_volume_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_brain_volume_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
