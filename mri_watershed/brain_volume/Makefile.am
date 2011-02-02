## 
## Makefile.am 
##

AM_CXXFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_brain_volume
mri_brain_volume_SOURCES=mri_brain_volume.cpp mrivariables.h TVector.h
mri_brain_volume_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_brain_volume_LDFLAGS=$(OS_LDFLAGS)

check_PROGRAMS = test
test_SOURCES=test.cpp TVector.h
test_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
test_LDFLAGS=$(OS_LDFLAGS)

TESTS = test

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
