## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_annotation2label
mri_annotation2label_SOURCES=mri_annotation2label.c
mri_annotation2label_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_annotation2label_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST=runit

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
