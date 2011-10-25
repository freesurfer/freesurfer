##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_ms_EM_with_atlas
mri_ms_EM_with_atlas_SOURCES=mri_ms_EM_with_atlas.c \
	PoissonSolver.c PoissonSolver.h 
mri_ms_EM_with_atlas_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_ms_EM_with_atlas_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
