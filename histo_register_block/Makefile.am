## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = histo_register_block
histo_register_block_SOURCES=histo_register_block.c
histo_register_block_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
histo_register_block_LDFLAGS=$(OS_LDFLAGS)

