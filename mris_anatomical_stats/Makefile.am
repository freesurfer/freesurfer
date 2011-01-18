## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_anatomical_stats
mris_anatomical_stats_SOURCES=mris_anatomical_stats.c
mris_anatomical_stats_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_anatomical_stats_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_anatomical_stats.help.xml

# trick to get test data into build dir
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mris_anatomical_stats/test1.expected"; then \
		cp -v $(top_srcdir)/mris_anatomical_stats/*.expected \
					$(top_builddir)/mris_anatomical_stats/; fi;

check_PROGRAMS=foo

TESTS=test_mris_anatomical_stats

EXTRA_DIST=test_mris_anatomical_stats test1.expected test2.expected \
test2_logfile.expected test3.expected test4.expected $(foo_DATA)

clean:
	rm -f mris_anatomical_stats *.o
	rm -Rf testdata2 *.observed *.table *.out

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
