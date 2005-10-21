## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mris_distance_to_label
mris_distance_to_label_SOURCES=mris_distance_to_label.cpp
mris_distance_to_label_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH) utilscpp/libutilscpp.a )
mris_distance_to_label_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
