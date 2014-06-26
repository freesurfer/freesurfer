##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -pg
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_dct_align_binary mri_dct_align
mri_dct_align_binary_SOURCES=mri_dct_align_binary.c 
mri_dct_align_binary_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_dct_align_binary_LDFLAGS=$(OS_LDFLAGS)

mri_dct_align_SOURCES=mri_dct_align.c 
mri_dct_align_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_dct_align_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
