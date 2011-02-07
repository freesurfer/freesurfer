## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include $(PETSC_CFLAGS) \
	-I$(top_srcdir)/tetgen/include
AM_LDFLAGS=$(PETSC_LIBS)

# if there is a need to activate petsc support, delete the #PETSC text

#PETSC if HAVE_PETSC_LIBS
#PETSC LIBS_PETSC=-lpetsc -lmpiuni
#PETSC else
LIBS_PETSC=
#PETSC endif

bin_PROGRAMS = mris_register
mris_register_SOURCES=mris_register.c
mris_register_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	-L$(top_builddir)/tetgen -ltetgen $(LIBS_PETSC)
mris_register_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_register.help.xml


#PETSC if HAVE_PETSC_LIBS
#PETSC install-exec-hook:
#PETSC # put a wrapper around the bin, used to setup Petsc enviro vars
#PETSC 	cp $(srcdir)/mris_register $(DESTDIR)$(bindir)/mris_register.bin
#PETSC 	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/mris_register
#PETSC 	echo "source \$$FREESURFER_HOME/bin/petsc_setup" \
#PETSC 		>> $(DESTDIR)$(bindir)/mris_register
#PETSC 	echo "mris_register.bin \$$argv" >> $(DESTDIR)$(bindir)/mris_register
#PETSC 	chmod a+x $(DESTDIR)$(bindir)/mris_register
#PETSC uninstall-hook:
#PETSC 	rm -f $(DESTDIR)$(bindir)/mris_register.bin
#PETSC else
#PETSC install-exec-hook:
#PETSC endif

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=make_mris_register_movie.tcl $(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_register.help.xml.h

