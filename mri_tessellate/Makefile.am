## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_tessellate ico_supersample
mri_tessellate_SOURCES=mri_tessellate.c
mri_tessellate_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_tessellate_LDFLAGS=$(OS_LDFLAGS)

ico_supersample_SOURCES=ico_supersample.c
ico_supersample_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
ico_supersample_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
