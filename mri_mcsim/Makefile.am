##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_mcsim
mri_mcsim_SOURCES=mri_mcsim.c
mri_mcsim_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_mcsim_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
