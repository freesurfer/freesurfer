##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_train_autoencoder mri_apply_autoencoder
mri_train_autoencoder_SOURCES=mri_train_autoencoder.c
mri_train_autoencoder_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_train_autoencoder_LDFLAGS=$(OS_LDFLAGS)

mri_apply_autoencoder_SOURCES=mri_apply_autoencoder.c
mri_apply_autoencoder_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_apply_autoencoder_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
