##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_compute_volume_fractions mri_compute_layer_fractions mris_compute_layer_intensities mri_compute_volume_intensities
mri_compute_volume_fractions_SOURCES=mri_compute_volume_fractions.c
mri_compute_volume_fractions_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_compute_volume_fractions_LDFLAGS=$(OS_LDFLAGS)

mri_compute_volume_intensities_SOURCES=mri_compute_volume_intensities.c
mri_compute_volume_intensities_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_compute_volume_intensities_LDFLAGS=$(OS_LDFLAGS)

mri_compute_layer_fractions_SOURCES=mri_compute_layer_fractions.c
mri_compute_layer_fractions_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_compute_layer_fractions_LDFLAGS=$(OS_LDFLAGS)

mris_compute_layer_intensities_SOURCES=mris_compute_layer_intensities.c
mris_compute_layer_intensities_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_compute_layer_intensities_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
