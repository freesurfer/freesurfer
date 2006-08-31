## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom $(GSL_CFLAGS)
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_glmfit
mri_glmfit_SOURCES=mri_glmfit.c
mri_glmfit_LDADD= $(top_srcdir)/fsgdf/libfsgdf.a \
$(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_glmfit_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
