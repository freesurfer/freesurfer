##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_rbvpvc
mri_rbvpvc_SOURCES=mri_rbvpvc.c
mri_rbvpvc_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_rbvpvc_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
