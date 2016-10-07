## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include $(PETSC_CFLAGS)
AM_LDFLAGS=$(PETSC_LIBS)

# if there is a need to activate petsc support, delete the #PETSC text

#PETSC if HAVE_PETSC_LIBS
#PETSC LIBS_PETSC=-lpetsc -lmpiuni
#PETSC else
LIBS_PETSC=
#PETSC endif

bin_PROGRAMS = mris_congeal
mris_congeal_SOURCES=mris_congeal.c
mris_congeal_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_congeal_LDFLAGS=$(OS_LDFLAGS)


#PETSC if HAVE_PETSC_LIBS
#PETSC install-exec-hook:
#PETSC # put a wrapper around the bin, used to setup Petsc enviro vars
#PETSC 	cp $(srcdir)/mris_congeal $(DESTDIR)$(bindir)/mris_congeal.bin
#PETSC 	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/mris_congeal
#PETSC 	echo "source \$$FREESURFER_HOME/bin/petsc_setup" \
#PETSC 		>> $(DESTDIR)$(bindir)/mris_congeal
#PETSC 	echo "mris_congeal.bin \$$argv" >> $(DESTDIR)$(bindir)/mris_congeal
#PETSC 	chmod a+x $(DESTDIR)$(bindir)/mris_congeal
#PETSC uninstall-hook:
#PETSC 	rm -f $(DESTDIR)$(bindir)/mris_congeal.bin
#PETSC else
#PETSC install-exec-hook:
#PETSC endif

EXTRA_DIST=render_surface.tcl make_movie.csh

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
