##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_concatenate_lta
mri_concatenate_lta_SOURCES=mri_concatenate_lta.c
mri_concatenate_lta_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_concatenate_lta_LDFLAGS=$(OS_LDFLAGS)

check_PROGRAMS=

TESTS=test_mri_concatenate_lta

EXTRA_DIST=test_mri_concatenate_lta nick13_to_nickbase.lta \
	reference-nickbase_to_nick13.lta

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
