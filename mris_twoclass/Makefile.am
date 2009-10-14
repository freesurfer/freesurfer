## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_twoclass
mris_twoclass_SOURCES=mris_twoclass.c
mris_twoclass_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_twoclass_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
