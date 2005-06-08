## 
# Makefile.am - specify in addition to the standard
##

AM_CFLAGS=$(GLUT_CFLAGS) -I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_probedicom
mri_probedicom_SOURCES=mri_probedicom.c
mri_probedicom_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(GLUT_LIBS)
mri_probedicom_LDFLAGS=

