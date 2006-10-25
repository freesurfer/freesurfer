## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_glm
mris_glm_SOURCES=mris_glm.c
mris_glm_LDADD= $(top_srcdir)/fsgdf/libfsgdf.a \
	$(addprefix $(top_srcdir)/, $(LIBS_MGH)) 
mris_glm_LDFLAGS=

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
