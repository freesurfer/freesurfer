##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_hausdorff_dist
mris_hausdorff_dist_SOURCES=mris_hausdorff_dist.c 
mris_hausdorff_dist_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_hausdorff_dist_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
