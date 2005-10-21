## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_info
mris_info_SOURCES=mris_info.cpp
mris_info_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_info_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
