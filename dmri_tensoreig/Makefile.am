
AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = dmri_tensoreig
dmri_tensoreig_SOURCES=dmri_tensoreig.c
dmri_tensoreig_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
dmri_tensoreig_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST=do.sh

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
