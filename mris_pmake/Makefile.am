##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_pmake
mris_pmake_SOURCES = \
	mris_pmake.cpp \
	abs2rel.c \
	c_label.cpp \
	c_SSocket.cpp \
	c_vertex.cpp \
	env.cpp \
	help.cpp \
	rel2abs.c \
	surface.cpp \
	asynch.cpp \
	c_SMessage.cpp \
	c_surface.cpp \
	dijkstra.cpp \
	general.cpp \
	scanopt.cpp
mris_pmake_LDADD = $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_pmake_LDFLAGS = $(OS_LDFLAGS)

EXTRA_DIST=asynch.h c_label.h c_SMessage.h c_SSocket.h c_surface.h \
	c_vertex.h dijkstra.h env.h general.h help.h pathconvert.h pstream.h \
	scanopt.h surface.h

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
