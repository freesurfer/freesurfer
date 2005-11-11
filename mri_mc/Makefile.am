AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_mc		
mri_mc_SOURCES=mri_mc.c
mri_mc_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
#mri_mc_LDFLAGS=$(OS_LDFLAGS)

check_PROGRAMS=myown build_test_volume

myown_SOURCES=myown.c
myown_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))
 
build_test_volume_SOURCES=build_test_volume.c
build_test_volume_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
