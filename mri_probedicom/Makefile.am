## 
# Makefile.am
##

AM_CFLAGS=$(GLUT_CFLAGS) -I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom
AM_LDFLAGS=

if HAVE_GLUT_LIBS
bin_PROGRAMS = mri_probedicom
mri_probedicom_SOURCES=mri_probedicom.c
mri_probedicom_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(GLUT_LIBS)
mri_probedicom_LDFLAGS=
else
bin_PROGRAMS=
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
