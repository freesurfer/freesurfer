##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_calc
mris_calc_SOURCES=mris_calc.c
mris_calc_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_calc_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
