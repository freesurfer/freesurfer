## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_hires_register \
	mri_nl_align_binary \
	mri_linear_align_binary \
	mri_compute_distances \
	mri_linear_align \
	mri_nl_align

mri_hires_register_SOURCES=mri_hires_register.c
mri_hires_register_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_hires_register_LDFLAGS=$(OS_LDFLAGS)

mri_nl_align_binary_SOURCES=mri_nl_align_binary.c
mri_nl_align_binary_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_nl_align_binary_LDFLAGS=$(OS_LDFLAGS)

mri_linear_align_binary_SOURCES=mri_linear_align_binary.c
mri_linear_align_binary_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_linear_align_binary_LDFLAGS=$(OS_LDFLAGS)

mri_compute_distances_SOURCES=mri_compute_distances.c
mri_compute_distances_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_compute_distances_LDFLAGS=$(OS_LDFLAGS)

mri_nl_align_SOURCES=mri_nl_align.c
mri_nl_align_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_nl_align_LDFLAGS=$(OS_LDFLAGS)

mri_linear_align_SOURCES=mri_linear_align.c
mri_linear_align_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_linear_align_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST = make_movie.tcl

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
