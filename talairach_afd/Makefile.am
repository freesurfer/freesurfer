##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = talairach_afd
talairach_afd_SOURCES=talairach_afd.c
talairach_afd_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
talairach_afd_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=talairach_afd.help.xml

# trick to get test data into build dir
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/talairach_afd/good_tal.txt"; then \
		cp -v $(top_srcdir)/talairach_afd/*tal.*  \
					$(top_builddir)/talairach_afd/; fi;
	if ! test -f "$(top_builddir)/distribution/fsafd/PonsCutDice.adf"; then \
		cp -v $(top_srcdir)/distribution/fsafd/*.adf  \
					$(top_builddir)/distribution/fsafd/; fi;

check_PROGRAMS=foo

TESTS=test_tal_afd

EXTRA_DIST=test_tal_afd good_tal.xfm good_tal.txt bad_tal.xfm bad_tal.txt \
	$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
