## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_cor2label
mri_cor2label_SOURCES=mri_cor2label.c
mri_cor2label_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_cor2label_LDFLAGS=$(OS_LDFLAGS)

# mri_vol2label is the better name
install-exec-hook:mri_cor2label
	cp $(top_builddir)/mri_cor2label/mri_cor2label \
		$(DESTDIR)$(bindir)/mri_vol2label

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
