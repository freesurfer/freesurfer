## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_segreg
mri_segreg_SOURCES=mri_segreg.c
mri_segreg_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_segreg_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
