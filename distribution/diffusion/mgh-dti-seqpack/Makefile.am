#
# Makefile.am
#

foodir=$(prefix)/diffusion/mgh-dti-seqpack
foo_DATA=$(wildcard *.gdt) README

EXTRA_DIST=$(foo_DATA)

install: 
	rm   -Rf $(DESTDIR)$(prefix)/diffusion/mgh-dti-seqpack/
	mkdir -p $(DESTDIR)$(prefix)/diffusion/mgh-dti-seqpack/
	cp $(foo_DATA) $(DESTDIR)$(prefix)/diffusion/mgh-dti-seqpack/

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
