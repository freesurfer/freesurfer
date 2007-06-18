## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include $(PETSC_CFLAGS) \
	-I$(top_srcdir)/tetgen/include
AM_LDFLAGS=$(PETSC_LIBS)

if HAVE_PETSC_LIBS
LIBS_PETSC=-lpetsc -lmpiuni
else
LIBS_PETSC=
endif

bin_PROGRAMS = mris_register
mris_register_SOURCES=mris_register.c
mris_register_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) \
	-L$(top_srcdir)/tetgen -ltetgen $(LIBS_PETSC)
mris_register_LDFLAGS=$(OS_LDFLAGS)

if HAVE_PETSC_LIBS
# put a wrapper around the bin, used to setup Petsc enviro vars
install-exec-hook:
	cp $(srcdir)/mris_register $(DESTDIR)$(bindir)/mris_register.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/mris_register
	echo "source \$$FREESURFER_HOME/bin/petsc_setup" \
		>> $(DESTDIR)$(bindir)/mris_register
	echo "mris_register.bin \$$argv" >> $(DESTDIR)$(bindir)/mris_register
	chmod a+x $(DESTDIR)$(bindir)/mris_register
endif

EXTRA_DIST=make_mris_register_movie.tcl

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
