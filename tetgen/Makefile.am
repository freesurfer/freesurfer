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

check_PROGRAMS = tetgen

TESTS = test_tetgen

EXTRA_DIST = test_tetgen example.poly \
test_example.1.face test_example.1.ele test_example.1.node \
README LICENSE README-MGH

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
