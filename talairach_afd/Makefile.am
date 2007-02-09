##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = talairach_afd
talairach_afd_SOURCES=talairach_afd.c
talairach_afd_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
talairach_afd_LDFLAGS=$(OS_LDFLAGS)

TESTS=test_tal_afd

EXTRA_DIST=test_tal_afd good_tal.xfm good_tal.txt bad_tal.xfm bad_tal.txt

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
