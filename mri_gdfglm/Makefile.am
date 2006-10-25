## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_gdfglm
mri_gdfglm_SOURCES=mri_gdfglm.c
mri_gdfglm_LDADD=$(top_srcdir)/fsgdf/libfsgdf.a \
	$(addprefix $(top_srcdir)/, $(LIBS_MGH)) 
mri_gdfglm_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
