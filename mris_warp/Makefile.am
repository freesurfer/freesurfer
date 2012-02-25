##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_warp
mris_warp_SOURCES=mris_warp.c
mris_warp_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_warp_LDFLAGS=$(OS_LDFLAGS)

#check_PROGRAMS=mris_warp

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
