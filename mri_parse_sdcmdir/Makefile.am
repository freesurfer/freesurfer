## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS=mri_parse_sdcmdir
mri_parse_sdcmdir_SOURCES=mri_parse_sdcmdir.c
mri_parse_sdcmdir_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_parse_sdcmdir_LDFLAGS=$(OS_LDFLAGS)

TESTS=test_mri_parse_sdcmdir

EXTRA_DIST=test_mri_parse_sdcmdir 556000-3-1.dcm test_dicomdir.sumfile 

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
