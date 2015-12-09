## 
# Makefile.am
##

AM_CFLAGS=-I$(top_srcdir)/glut $(GL_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom
AM_LDFLAGS=

if HAVE_OPENGL_LIBS
bin_PROGRAMS = mri_probedicom
mri_probedicom_SOURCES=mri_probedicom.c
mri_probedicom_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(LIBS_OPENGL) $(top_builddir)/glut/libglut.a $(X_LIBS)
mri_probedicom_LDFLAGS=
else
bin_PROGRAMS=
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
