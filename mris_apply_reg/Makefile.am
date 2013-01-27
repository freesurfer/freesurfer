##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_apply_reg
mris_apply_reg_SOURCES=mris_apply_reg.c mris_apply_reg.help.xml.h
mris_apply_reg_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_apply_reg_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_apply_reg.help.xml

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_apply_reg.help.xml.h

