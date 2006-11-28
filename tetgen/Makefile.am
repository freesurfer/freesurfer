#
# Makefile.am for TetGen
#

# build library...
noinst_LIBRARIES= libtetgen.a
libtetgen_a_SOURCES = predicates.cxx tetgen.cxx
libtetgen_a_CXXFLAGS = -DSELF_CHECK -DTETLIBRARY

# and binary:
bin_PROGRAMS = tetgen
tetgen_SOURCES = predicates.cxx tetgen.cxx
tetgen_CXXFLAGS = -DSELF_CHECK

TESTS = test_tetgen

EXTRA_DIST = test_tetgen example.poly \
test_example.1.face test_example.1.ele test_example.1.node

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
