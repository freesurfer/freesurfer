##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_calc
mris_calc_SOURCES=mris_calc.c
mris_calc_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_calc_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_calc.help.xml

EXTRA_DIST=$(foo_DATA)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
