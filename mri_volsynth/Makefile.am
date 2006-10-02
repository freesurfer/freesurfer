## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_volsynth
mri_volsynth_SOURCES=mri_volsynth.c
mri_volsynth_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_volsynth_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
