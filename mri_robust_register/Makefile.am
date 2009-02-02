##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_robust_register
mri_robust_register_SOURCES=mri_robust_register.cpp Regression.cpp Registration.cpp
mri_robust_register_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_robust_register_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST=Regression.h Registration.h Quaternion.h

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
