##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_aseg_distance
mris_aseg_distance_SOURCES=mris_aseg_distance.c
mris_aseg_distance_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_aseg_distance_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
