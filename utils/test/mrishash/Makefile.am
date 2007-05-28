## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

check_PROGRAMS = mrishash_demo_100_find_coverage \
	mrishash_demo_200_mht_hatch  \
	mrishash_test_100_find_tests \
	mrishash_test_200_intersect

TESTS=mrishash_test_100_find_tests mrishash_test_200_intersect

#------------- exercise ----------------

mrishash_demo_100_find_coverage_SOURCES=mrishash_demo_100_find_coverage.c
mrishash_demo_100_find_coverage_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mrishash_demo_100_find_coverage_LDFLAGS= $(OS_LDFLAGS)

mrishash_demo_200_mht_hatch_SOURCES=mrishash_demo_200_mht_hatch.c
mrishash_demo_200_mht_hatch_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mrishash_demo_200_mht_hatch_LDFLAGS= $(OS_LDFLAGS)

##------------- test ----------------

mrishash_test_100_find_tests_SOURCES=mrishash_test_100_find_tests.c
mrishash_test_100_find_tests_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mrishash_test_100_find_tests_LDFLAGS= $(OS_LDFLAGS)

mrishash_test_200_intersect_SOURCES=mrishash_test_200_intersect.c
mrishash_test_200_intersect_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mrishash_test_200_intersect_LDFLAGS= $(OS_LDFLAGS)

EXTRA_DIST=

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	rm -f *_mri.mgz *_log.txt *.tri *.o
