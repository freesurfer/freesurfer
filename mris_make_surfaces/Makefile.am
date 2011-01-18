## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_make_surfaces mris_refine_surfaces mris_mef_surfaces mris_exvivo_surfaces
mris_make_surfaces_SOURCES=mris_make_surfaces.c
mris_make_surfaces_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_make_surfaces_LDFLAGS= $(OS_LDFLAGS)

mris_refine_surfaces_SOURCES=mris_refine_surfaces.c
mris_refine_surfaces_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_refine_surfaces_LDFLAGS= $(OS_LDFLAGS)

mris_mef_surfaces_SOURCES=mris_mef_surfaces.c
mris_mef_surfaces_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_mef_surfaces_LDFLAGS= $(OS_LDFLAGS)

mris_exvivo_surfaces_SOURCES=mris_exvivo_surfaces.c
mris_exvivo_surfaces_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_exvivo_surfaces_LDFLAGS= $(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_make_surfaces.help.xml

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
