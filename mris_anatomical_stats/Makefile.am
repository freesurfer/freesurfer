## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_anatomical_stats
mris_anatomical_stats_SOURCES=mris_anatomical_stats.c
mris_anatomical_stats_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_anatomical_stats_LDFLAGS=$(OS_LDFLAGS)

TESTS=test_mris_anatomical_stats

EXTRA_DIST=test_mris_anatomical_stats

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
