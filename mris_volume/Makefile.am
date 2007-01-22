## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mris_volume
mris_volume_SOURCES=mris_volume.c
mris_volume_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_volume_LDFLAGS=$(OS_LDFLAGS)


bin_PROGRAMS = mris_wm_volume
mris_wm_volume_SOURCES=mris_wm_volume.c
mris_wm_volume_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_wm_volume_LDFLAGS=$(OS_LDFLAGS)




# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
