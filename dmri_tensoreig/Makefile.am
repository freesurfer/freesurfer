
AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

if NMR_INSTALL
bin_PROGRAMS = dmri_tensoreig
dmri_tensoreig_SOURCES=dmri_tensoreig.c
dmri_tensoreig_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
dmri_tensoreig_LDFLAGS=$(OS_LDFLAGS)
else
bin_PROGRAMS=
endif

EXCLUDE_FILES=dmri_tensoreig.c do.sh
include $(top_srcdir)/Makefile.extra
