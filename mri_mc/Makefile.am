AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_mc mri_pretess mris_extract_main_component

mri_mc_SOURCES=mri_mc.c
mri_mc_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))

foodir=$(prefix)/docs/xml
foo_DATA=mri_pretess.help.xml mris_extract_main_component.help.xml

# trick to get test data into build dir
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_mc/test_volume_reference.mgz"; then \
		cp -v $(top_srcdir)/mri_mc/test_volume_reference.mgz \
					$(top_builddir)/mri_mc/; fi;

check_PROGRAMS=foo build_test_volume
#build_MC_table

TESTS=compare_test_volumes

build_test_volume_SOURCES=build_test_volume.c
build_test_volume_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH))

mri_pretess_SOURCES=mri_pretess.c
mri_pretess_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH))

mris_extract_main_component_SOURCES=mris_extract_main_component.c
mris_extract_main_component_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH))

#build_MC_table_SOURCES=build_MC_table.c
#build_MC_table_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH))

EXTRA_DIST=compare_test_volumes test_volume_reference.mgz $(foo_DATA) \
	build_MC_table.c mc_table.gif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
