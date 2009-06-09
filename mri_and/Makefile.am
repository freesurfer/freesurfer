##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_and mri_or

mri_and_SOURCES=mri_and.c
mri_and_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_and_LDFLAGS=$(OS_LDFLAGS)

mri_or_SOURCES=mri_or.c
mri_or_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_or_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
