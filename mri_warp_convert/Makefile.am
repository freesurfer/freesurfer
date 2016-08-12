##
## Makefile.am
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

BUILT_SOURCES=mri_warp_convert.help.xml.h mri_warp_convert.help.xml.html

foodir=$(prefix)/docs/xml
foo_DATA=mri_warp_convert.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mri_warp_convert.help.xml.html

bin_PROGRAMS = mri_warp_convert
mri_warp_convert_SOURCES=mri_warp_convert.cpp
mri_warp_convert_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_warp_convert_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	rm -f $(BUILT_SOURCES)
