## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_train
mri_train_SOURCES=mri_train.c
mri_train_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_train_LDFLAGS=$(OS_LDFLAGS)

