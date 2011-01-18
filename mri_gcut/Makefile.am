AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_gcut
mri_gcut_SOURCES=mri_gcut.cpp
mri_gcut_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_gcut_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_gcut.help.xml

EXTRA_DIST = gcut.h graphcut.h graphcut.cpp pre_pro.cpp $(foo_DATA)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
