## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_register_to_volume mris_register_to_label
mris_register_to_volume_SOURCES=mris_register_to_volume.c
mris_register_to_volume_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_register_to_volume_LDFLAGS=$(OS_LDFLAGS)

mris_register_to_label_SOURCES=mris_register_to_label.c
mris_register_to_label_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_register_to_label_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
