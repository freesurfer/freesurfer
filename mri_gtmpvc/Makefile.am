##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_gtmpvc
mri_gtmpvc_SOURCES=mri_gtmpvc.c
mri_gtmpvc_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_gtmpvc_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
