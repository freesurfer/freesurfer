## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_tessellate ico_supersample
mri_tessellate_SOURCES=mri_tessellate.c
mri_tessellate_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_tessellate_LDFLAGS=$(OS_LDFLAGS)

ico_supersample_SOURCES=ico_supersample.c trilib.h systemid.h
ico_supersample_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
ico_supersample_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_tessellate.help.xml

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
