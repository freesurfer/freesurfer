## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_ca_train
mris_ca_train_SOURCES=mris_ca_train.c
mris_ca_train_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_ca_train_LDFLAGS=$(OS_LDFLAGS)

