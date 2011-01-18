#
# Makefile.am for TetGen
#

# need -O0 flag, as 'tetgen example.poly' will coredump on 32bit platfomrs 
# (but not 64bit!) if compiled with optimization enabled

# build library...
noinst_LIBRARIES= libtetgen.a
libtetgen_a_SOURCES = predicates.cxx tetgen.cxx
libtetgen_a_CXXFLAGS = -DSELF_CHECK -DTETLIBRARY -O0

# and binary:
bin_PROGRAMS = tetgen
tetgen_SOURCES = predicates.cxx tetgen.cxx
tetgen_CXXFLAGS = -DSELF_CHECK -O0

# trick to get test data into build dir
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/tetgen/example.poly"; then \
		cp -v $(top_srcdir)/tetgen/example.poly $(top_srcdir)/tetgen/*.1.*  \
					$(top_builddir)/tetgen/; fi;

check_PROGRAMS = foo tetgen

TESTS = test_tetgen

EXTRA_DIST = tetgen.h test_tetgen example.poly \
	x86_64-example.1.face x86_64-example.1.ele x86_64-example.1.node \
	README LICENSE README-MGH

DISTCLEANFILES=example.1.face example.1.node example.1.ele

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
