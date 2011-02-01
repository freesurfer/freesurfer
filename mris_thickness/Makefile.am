## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_thickness mris_intensity_profile
mris_thickness_SOURCES=mris_thickness.c
mris_thickness_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_thickness_LDFLAGS=$(OS_LDFLAGS)

mris_intensity_profile_SOURCES=mris_intensity_profile.c
mris_intensity_profile_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_intensity_profile_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST = mris_cluster_profiles.c

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
