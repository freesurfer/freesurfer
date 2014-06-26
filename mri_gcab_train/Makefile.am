##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_gcab_train
mri_gcab_train_SOURCES=mri_gcab_train.c
mri_gcab_train_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_gcab_train_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
