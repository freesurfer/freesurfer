AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_mc mri_pretess mris_extract_main_component

mri_mc_SOURCES=mri_mc.c
mri_mc_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
#mri_mc_LDFLAGS=$(OS_LDFLAGS)

check_PROGRAMS=build_test_volume mri_pretess mris_extract_main_component
#build_MC_table

build_test_volume_SOURCES=build_test_volume.c
build_test_volume_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

mri_pretess_SOURCES=mri_pretess.c
mri_pretess_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

mris_extract_main_component_SOURCES=mris_extract_main_component.c
mris_extract_main_component_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

#build_MC_table_SOURCES=build_MC_table.c
#build_MC_table_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
