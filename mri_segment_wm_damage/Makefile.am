##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_segment_wm_damage
mri_segment_wm_damage_SOURCES=mri_segment_wm_damage.c
mri_segment_wm_damage_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_segment_wm_damage_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
