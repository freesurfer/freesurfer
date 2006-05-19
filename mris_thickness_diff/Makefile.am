##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_thickness_diff
mris_thickness_diff_SOURCES=mris_thickness_diff.c
mris_thickness_diff_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_thickness_diff_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
