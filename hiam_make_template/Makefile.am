##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = hiam_make_template
hiam_make_template_SOURCES=hiam_make_template.c
hiam_make_template_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
hiam_make_template_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

