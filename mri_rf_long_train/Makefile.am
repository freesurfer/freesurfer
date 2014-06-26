##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_rf_long_train
mri_rf_long_train_SOURCES=mri_rf_long_train.c
mri_rf_long_train_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_rf_long_train_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
